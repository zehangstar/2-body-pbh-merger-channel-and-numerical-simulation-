"""Binary-single 的单时间片、单壳层批量轨道演化。

轨道方程消费 ``HaloShellState`` 中的环境密度、速度与硬双星边界；
``shell_density_msun_pc3`` 始终是未缩放的 NFW 局部密度。
这里仅返回原始 Monte Carlo 事件，不进行人口重加权或 HMF 积分。

论文原文印刷 Eq. (24) 对应 mu=2；Table III 诊断更接近 mu=0.5。
本项目以 Prada12-HMF、mu=0.5、壳层中点速度、ksi=1 为默认配置，
任何表格匹配分支必须在调用时显式指定，不能混称论文原始模型。
"""

from dataclasses import dataclass, replace
import math

import numpy as np
from scipy.integrate import RK45

import haloconcentration as hc
import halostructure as hs


G_SI = 6.67430e-11
C_SI = 299792458.0
SOLAR_MASS_KG = 1.98847e30
PARSEC_M = 3.085677581491367e16
MYR_S = 365.25 * 86400.0 * 1.0e6

# Sesana, Haardt & Madau (2006), ApJ 651, 392, Eq. (18), Table 3,
# https://arxiv.org/pdf/astro-ph/0604299 . q=M2/M1=1, e=0.15...0.90.
# a0 列的单位是 a_h；K=A(1+a/(a0*a_h))**gamma+B。
SESANA_EQUAL_MASS_E = np.array([0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90])
SESANA_EQUAL_MASS_A = np.array([0.0, 0.037, 0.075, 0.105, 0.121, 0.134, 0.082])
SESANA_EQUAL_MASS_A0_OVER_AH = np.array([1.0, 0.339, 0.151, 0.088, 0.090, 0.064, 0.085])
SESANA_EQUAL_MASS_GAMMA = np.array([0.0, -3.335, -1.548, -0.893, -0.895, -0.544, -0.663])
SESANA_EQUAL_MASS_B = np.array([0.0, -0.012, -0.008, -0.005, -0.008, -0.006, -0.004])


@dataclass(frozen=True)
class BinarySingleConfig:
    """后续轨道计算可复用的物理与时间步配置。

    ``ksi`` 定义 ``rho_env = ksi * rho_NFW(midpoint)``，不改变
    halo 质量、浓度或 NFW 剖面。``velocity_mu`` 定义
    ``v_env = sqrt(mu * G * M(<r) / r)``；``velocity_radius_strategy``
    独立选择该速度在壳中点或壳外边界评价。``ksi`` 与 PBH 的
    单体/双星人口比例是不同参数，不能互相代用。
    """

    halo_model: str = hc.DEFAULT_HALO_MODEL
    primary_mass_msun: float = 30.0
    secondary_mass_msun: float = 30.0
    pbh_fraction_total: float = 1.0
    fraction_single: float = 0.5
    fraction_binary: float = 0.5
    velocity_mu: float = 0.5
    velocity_radius_strategy: str = "midpoint"
    ksi: float = 1.0
    global_timestep_myr: float = 200.0
    local_timestep_myr: float = 2.0
    start_redshift: float = 12.0
    sample_count_per_shell_and_step: int = 2_000_000
    integration_method: str = "adaptive_log_a"
    eccentricity_growth_model: str = "sesana_endpoint"
    integration_rtol: float = 1.0e-6
    integration_time_atol_myr: float = 1.0e-8
    integration_log_j2_atol: float = 1.0e-9
    integration_max_log_a_step: float = 0.25
    integration_max_attempts: int = 10000
    merger_radius_gm_c2: float = 6.0

    def __post_init__(self):
        for name in (
            "primary_mass_msun",
            "secondary_mass_msun",
            "velocity_mu",
            "global_timestep_myr",
            "local_timestep_myr",
            "integration_rtol",
            "integration_time_atol_myr",
            "integration_log_j2_atol",
            "integration_max_log_a_step",
            "merger_radius_gm_c2",
        ):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} 必须为有限正数。")
        if not math.isfinite(self.ksi) or self.ksi < 0.0:
            raise ValueError("ksi 必须为有限非负数。")
        if not 0.0 <= self.pbh_fraction_total <= 1.0:
            raise ValueError("pbh_fraction_total 必须在 [0, 1] 内。")
        if not 0.0 <= self.fraction_single <= 1.0:
            raise ValueError("fraction_single 必须在 [0, 1] 内。")
        if not 0.0 <= self.fraction_binary <= 1.0:
            raise ValueError("fraction_binary 必须在 [0, 1] 内。")
        if not math.isclose(
            self.fraction_single + self.fraction_binary, 1.0,
            rel_tol=0.0, abs_tol=1e-12,
        ):
            raise ValueError("fraction_single 与 fraction_binary 之和必须为 1。")
        if not math.isfinite(self.start_redshift) or self.start_redshift < 0.0:
            raise ValueError("start_redshift 必须为有限非负数。")
        if self.sample_count_per_shell_and_step <= 0:
            raise ValueError("sample_count_per_shell_and_step 必须为正整数。")
        if self.integration_method not in ("euler", "adaptive_log_a"):
            raise ValueError("integration_method 必须是 euler 或 adaptive_log_a。")
        if self.velocity_radius_strategy not in hs.VELOCITY_RADIUS_STRATEGIES:
            choices = ", ".join(hs.VELOCITY_RADIUS_STRATEGIES)
            raise ValueError(f"velocity_radius_strategy 必须是 {choices} 之一。")
        if self.eccentricity_growth_model not in ("sesana_endpoint", "sesana_taper", "zero"):
            raise ValueError("未知 eccentricity_growth_model。")
        if not isinstance(self.integration_max_attempts, int) or self.integration_max_attempts <= 0:
            raise ValueError("integration_max_attempts 必须为正整数。")
        steps = self.global_timestep_myr / self.local_timestep_myr
        if self.integration_method == "euler" and not math.isclose(steps, round(steps), rel_tol=0.0, abs_tol=1e-10):
            raise ValueError("全局时间步必须是局部时间步的整数倍。")

    @property
    def local_steps_per_global_step(self):
        """仅用于 Euler 分支；论文基准为 200/2=100 个局部步。"""
        return round(self.global_timestep_myr / self.local_timestep_myr)


def build_shell_environment(present_day_halo_mass_msun, z, config=None):
    """由固定 ``M0`` 轨迹构造某个红移的完整壳层环境。

    返回 ``halostructure.HaloShellState``，数组下标从 0 开始，对应论文
    球壳编号 1、2、...。壳层数沿轨迹固定，边界随瞬时 ``R_vir(z)`` 更新。
    本函数不抽样双星，也不计算轨道微分方程或事件数。
    """
    if config is None:
        config = BinarySingleConfig()
    return hs.build_halo_shell_state(
        present_day_halo_mass_msun,
        z,
        primary_pbh_mass_msun=config.primary_mass_msun,
        model=config.halo_model,
        velocity_mu=config.velocity_mu,
        velocity_radius_strategy=config.velocity_radius_strategy,
        ksi=config.ksi,
    )


@dataclass(frozen=True)
class BinaryOrbitState:
    """单壳层中一批双星的当前轨道，数组单位分别为 pc 和无量纲。"""

    semi_major_axis_pc: np.ndarray
    eccentricity: np.ndarray
    active_mask: np.ndarray


@dataclass(frozen=True)
class OrbitalDerivatives:
    """Eqs. (19)、(25) 的四个贡献，单位为 pc/Myr 与 1/Myr。"""

    three_body_da_pc_myr: np.ndarray
    gravitational_wave_da_pc_myr: np.ndarray
    three_body_de_per_myr: np.ndarray
    gravitational_wave_de_per_myr: np.ndarray

    @property
    def da_pc_myr(self):
        return self.three_body_da_pc_myr + self.gravitational_wave_da_pc_myr

    @property
    def de_per_myr(self):
        return self.three_body_de_per_myr + self.gravitational_wave_de_per_myr


@dataclass(frozen=True)
class BinaryEvolutionEvents:
    """一个冻结环境时间片的原始事件；没有物理双星数重加权。

    adaptive_log_a 中 merger_step 是该轨道接受的积分步编号，
    completed_local_steps 是批内最大接受步数；二者不再代表 2 Myr 网格。
    已并合轨道 final_state.a=0 是终止标记，实际终点由配置的有限半径定义。
    """

    final_state: BinaryOrbitState
    hard_mask: np.ndarray
    merged_mask: np.ndarray
    invalid_mask: np.ndarray
    k_above_calibration_mask: np.ndarray
    merger_step: np.ndarray
    merger_time_myr: np.ndarray
    completed_local_steps: int

    @property
    def sample_count(self):
        return self.hard_mask.size

    @property
    def hard_count(self):
        return int(np.count_nonzero(self.hard_mask))

    @property
    def merger_count(self):
        return int(np.count_nonzero(self.merged_mask))

    @property
    def invalid_count(self):
        return int(np.count_nonzero(self.invalid_mask))


def sesana_equal_mass_k(semi_major_axis_pc, eccentricity, hard_semimajor_axis_pc):
    """Sesana et al. (2006) Eq. (18), Table 3 的等质量 ``K`` 拟合。

    原始散射实验是大质量黑洞双星与轻得多的场星，背景为固定、各向同性
    Maxwellian 分布；移用于 PBH 与同质量第三体是目标 PBH 论文的外推。
    表 3 只列 e=0.15..0.9 的离散初始偏心率。这里按当前 e 对 *K 值*
    分段线性插值，e=0 的 K=0；e>0.9 固定采用 0.9 的端点值。
    这两项插值/端点规则是本项目数值约定，不是 Sesana 原文公式。
    """
    a, e, ah = np.broadcast_arrays(
        np.asarray(semi_major_axis_pc, dtype=float),
        np.asarray(eccentricity, dtype=float),
        np.asarray(hard_semimajor_axis_pc, dtype=float),
    )
    if np.any(~np.isfinite(a)) or np.any(a <= 0.0):
        raise ValueError("半长轴必须为有限正数。")
    if np.any(~np.isfinite(e)) or np.any((e < 0.0) | (e >= 1.0)):
        raise ValueError("偏心率必须满足 0 <= e < 1。")
    if np.any(~np.isfinite(ah)) or np.any(ah <= 0.0):
        raise ValueError("硬双星半长轴必须为有限正数。")

    e_fit = np.minimum(e, SESANA_EQUAL_MASS_E[-1])
    upper = np.clip(
        np.searchsorted(SESANA_EQUAL_MASS_E, e_fit, side="right"),
        1, len(SESANA_EQUAL_MASS_E) - 1,
    )
    lower = upper - 1
    x = a / ah

    def fitted_at(index):
        return (
            SESANA_EQUAL_MASS_A[index]
            * (1.0 + x / SESANA_EQUAL_MASS_A0_OVER_AH[index])
            ** SESANA_EQUAL_MASS_GAMMA[index]
            + SESANA_EQUAL_MASS_B[index]
        )

    weight = (
        (e_fit - SESANA_EQUAL_MASS_E[lower])
        / (SESANA_EQUAL_MASS_E[upper] - SESANA_EQUAL_MASS_E[lower])
    )
    return fitted_at(lower) * (1.0 - weight) + fitted_at(upper) * weight


def hardening_coefficient(semi_major_axis_pc, hard_semimajor_axis_pc):
    """目标 PBH 论文 Eq. (22)，即 Sesana 等质量 H 的拟合。"""
    a = np.asarray(semi_major_axis_pc, dtype=float)
    ah = np.asarray(hard_semimajor_axis_pc, dtype=float)
    return 14.55 * (1.0 + 0.287 * a / ah) ** -0.95


def eccentricity_growth_coefficient(a, e, ah, model="sesana_endpoint"):
    """主设置为原有端点 K；taper 与 zero 仅用于明确标注的敏感性诊断。

    taper 在 e>0.9 乘以 (1-e²)/(1-0.9²)，不是新的 Sesana 拟合。
    """
    if model == "zero":
        return np.zeros_like(np.asarray(a, dtype=float))
    k = sesana_equal_mass_k(a, e, ah)
    if model == "sesana_taper":
        k = k * np.minimum(1.0, (1.0 - np.asarray(e)**2) / 0.19)
    elif model != "sesana_endpoint":
        raise ValueError("未知 eccentricity_growth_model。")
    return k


def peters_f(eccentricity):
    """目标 PBH 论文 Eq. (20)。"""
    e = np.asarray(eccentricity, dtype=float)
    return (1.0 + 73.0 * e**2 / 24.0 + 37.0 * e**4 / 96.0) / (1.0 - e**2) ** 3.5


def peters_d(eccentricity):
    """目标 PBH 论文 Eq. (26)。"""
    e = np.asarray(eccentricity, dtype=float)
    return (e + 121.0 * e**3 / 304.0) / (1.0 - e**2) ** 2.5


def orbital_derivatives(
    semi_major_axis_pc,
    eccentricity,
    environment_density_msun_pc3,
    velocity_dispersion_km_s,
    hard_semimajor_axis_pc,
    primary_mass_msun=30.0,
    secondary_mass_msun=30.0,
    eccentricity_growth_model="sesana_endpoint",
):
    """向量化计算目标 PBH 论文 Eqs. (19)、(25) 的分项右端。"""
    if primary_mass_msun != secondary_mass_msun:
        raise ValueError("当前 Sesana K 只实现 q=1；非等质量双星尚不可演化。")
    if primary_mass_msun <= 0.0 or not math.isfinite(primary_mass_msun):
        raise ValueError("PBH 质量必须为有限正数。")
    a, e = np.broadcast_arrays(
        np.asarray(semi_major_axis_pc, dtype=float),
        np.asarray(eccentricity, dtype=float),
    )
    rho = float(environment_density_msun_pc3)
    v = float(velocity_dispersion_km_s)
    ah = float(hard_semimajor_axis_pc)
    if not math.isfinite(rho) or rho < 0.0:
        raise ValueError("环境密度必须为有限非负数。")
    if not math.isfinite(v) or v <= 0.0:
        raise ValueError("速度弥散必须为有限正数。")
    k = eccentricity_growth_coefficient(a, e, ah, eccentricity_growth_model)
    h = hardening_coefficient(a, ah)

    a_m = a * PARSEC_M
    rho_kg_m3 = rho * SOLAR_MASS_KG / PARSEC_M**3
    velocity_m_s = v * 1000.0
    m1 = primary_mass_msun * SOLAR_MASS_KG
    m2 = secondary_mass_msun * SOLAR_MASS_KG
    gw_factor = G_SI**3 * (m1 + m2) * m1 * m2 / C_SI**5

    da_three_body = (
        -G_SI * h * rho_kg_m3 / velocity_m_s * a_m**2
        * MYR_S / PARSEC_M
    )
    da_gw = -64.0 / 5.0 * gw_factor / a_m**3 * peters_f(e) * MYR_S / PARSEC_M
    de_three_body = G_SI * h * k * rho_kg_m3 / velocity_m_s * a_m * MYR_S
    de_gw = -304.0 / 15.0 * gw_factor / a_m**4 * peters_d(e) * MYR_S
    return OrbitalDerivatives(da_three_body, da_gw, de_three_body, de_gw)


def evolve_binary_batch(
    initial_semi_major_axis_pc,
    initial_eccentricity,
    shell_state,
    shell_index,
    config=None,
):
    """在一个冻结的壳层环境中推进一批初态，返回原始并合事件。

    默认 adaptive_log_a 按容差控制独立轨道步长，见下方变量变换积分器。
    integration_method='euler' 显式选择原有论文步长诊断分支：
    显式 Euler 以 ``config.local_timestep_myr`` 推进一个
    ``config.global_timestep_myr`` 时间片。若 Euler 使 a<=0，记录为
    并合，并用该步线性插值估计事件时间；论文没有指定此越界规则。
    若 a 仍为正而 e 越出 [0,1) 或右端非有限，记录为数值失效，
    不把它静默计为并合。这里不继承样本到下一个全局时间片。
    """
    if config is None:
        config = BinarySingleConfig()
    a = np.array(initial_semi_major_axis_pc, dtype=float, copy=True)
    e = np.array(initial_eccentricity, dtype=float, copy=True)
    if a.ndim != 1 or e.shape != a.shape:
        raise ValueError("初始 a、e 必须是等长一维数组。")
    if np.any(~np.isfinite(a)) or np.any(a <= 0.0):
        raise ValueError("初始半长轴必须为有限正数。")
    if np.any(~np.isfinite(e)) or np.any((e < 0.0) | (e >= 1.0)):
        raise ValueError("初始偏心率必须满足 0 <= e < 1。")
    shell_index = int(shell_index)
    if shell_index < 0 or shell_index >= shell_state.shell_mass_msun.size:
        raise IndexError("shell_index 超出壳层范围。")
    if config.integration_method == "adaptive_log_a":
        return _evolve_adaptive_log_a(a, e, shell_state, shell_index, config)
    rho = float(shell_state.environment_density_msun_pc3[shell_index])
    v = float(shell_state.velocity_dispersion_km_s[shell_index])
    ah = float(shell_state.hard_semimajor_axis_pc[shell_index])
    hard = a <= ah
    active = hard.copy()
    merged = np.zeros(a.size, dtype=bool)
    invalid = np.zeros(a.size, dtype=bool)
    k_above_calibration = np.zeros(a.size, dtype=bool)
    merger_step = np.full(a.size, -1, dtype=int)
    merger_time_myr = np.full(a.size, np.nan)
    completed_steps = 0
    dt = config.local_timestep_myr

    for step in range(config.local_steps_per_global_step):
        active_index = np.flatnonzero(active)
        if active_index.size == 0:
            break
        k_above_calibration[active_index] |= e[active_index] > SESANA_EQUAL_MASS_E[-1]
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            rhs = orbital_derivatives(
                a[active_index], e[active_index], rho, v, ah,
                config.primary_mass_msun, config.secondary_mass_msun,
                config.eccentricity_growth_model,
            )
            next_a = a[active_index] + dt * rhs.da_pc_myr
            next_e = e[active_index] + dt * rhs.de_per_myr

        finite_rhs = np.isfinite(rhs.da_pc_myr) & np.isfinite(rhs.de_per_myr)
        newly_merged = finite_rhs & np.isfinite(next_a) & (next_a <= 0.0)
        newly_invalid = (
            ~finite_rhs | ~np.isfinite(next_a) | ~np.isfinite(next_e)
            | ((next_a > 0.0) & ((next_e < 0.0) | (next_e >= 1.0)))
        ) & ~newly_merged
        still_active = ~(newly_merged | newly_invalid)

        merged_index = active_index[newly_merged]
        if merged_index.size:
            fraction = a[merged_index] / (-dt * rhs.da_pc_myr[newly_merged])
            merger_time_myr[merged_index] = (step + fraction) * dt
            merger_step[merged_index] = step + 1
            a[merged_index] = 0.0
            merged[merged_index] = True
        invalid[active_index[newly_invalid]] = True
        survivor_index = active_index[still_active]
        a[survivor_index] = next_a[still_active]
        e[survivor_index] = next_e[still_active]
        active[active_index[newly_merged | newly_invalid]] = False
        completed_steps = step + 1

    return BinaryEvolutionEvents(
        final_state=BinaryOrbitState(a, e, active),
        hard_mask=hard,
        merged_mask=merged,
        invalid_mask=invalid,
        k_above_calibration_mask=k_above_calibration,
        merger_step=merger_step,
        merger_time_myr=merger_time_myr,
        completed_local_steps=completed_steps,
    )


def evolve_binary_population_batch(a, e, shell_state, shell_index, config=None):
    """推进 cohort 存活者：硬双星用完整方程，软双星仍受 GW 作用。

    返回的 hard_mask 是该片开始时的硬筛选；merged_mask 可以包含软双星
    的纯 GW 并合。软双星不加入环境散射、离解或反冲模型。调用者应在下片
    传入存活者的 final_state，不能重新抽取这些个体的初始轨道。
    """
    if config is None:
        config = BinarySingleConfig()
    events = evolve_binary_batch(a, e, shell_state, shell_index, config)
    soft = np.flatnonzero(~events.hard_mask)
    if soft.size:
        a_soft = np.asarray(a)[soft]
        isolated = replace(
            shell_state,
            environment_density_msun_pc3=np.zeros_like(shell_state.shell_mass_msun),
            hard_semimajor_axis_pc=np.full_like(shell_state.shell_mass_msun, np.max(a_soft)),
        )
        gw = evolve_binary_batch(
            a_soft, np.asarray(e)[soft], isolated, shell_index,
            replace(config, eccentricity_growth_model='zero'),
        )
        events.final_state.semi_major_axis_pc[soft] = gw.final_state.semi_major_axis_pc
        events.final_state.eccentricity[soft] = gw.final_state.eccentricity
        events.final_state.active_mask[soft] = ~gw.merged_mask & ~gw.invalid_mask
        events.merged_mask[soft] = gw.merged_mask
        events.invalid_mask[soft] = gw.invalid_mask
        events.merger_time_myr[soft] = gw.merger_time_myr
        events.merger_step[soft] = gw.merger_step
        events = replace(events, completed_local_steps=max(events.completed_local_steps, gw.completed_local_steps))
    return events


def log_a_orbital_rhs(u, state, initial_a_pc, rho, velocity, ah, config):
    """Eqs.19/25 的等价变量变换，返回 d(t, ln(j²))/du。

    u=ln(a_initial/a) 单调递增；j²=1-e² 单独保存，避免高 e 时
    相减丢失精度。GW 主导末期 du/dt 发散，但 dt/du 趋于零，
    dln(j²)/du 有限。此函数也用于独立 scipy DOP853 对照。
    非物理试探级返回 NaN，外层自适应积分器拒绝该步，不裁剪已接受轨道。
    """
    state = np.asarray(state, dtype=float)
    log_j2 = state[..., 1]
    valid = np.isfinite(log_j2) & (log_j2 <= 0.0) & (log_j2 > -700.0)
    safe_log_j2 = np.clip(np.nan_to_num(log_j2), -700.0, 0.0)
    e2 = -np.expm1(safe_log_j2)
    e = np.sqrt(e2)
    log_a = np.log(initial_a_pc) - u
    a = np.exp(log_a)
    k = eccentricity_growth_coefficient(a, np.minimum(e, np.nextafter(1.0, 0.0)), ah, config.eccentricity_growth_model)
    h = hardening_coefficient(a, ah)
    m1 = config.primary_mass_msun * SOLAR_MASS_KG
    m2 = config.secondary_mass_msun * SOLAR_MASS_KG
    beta = (64.0 / 5.0 * G_SI**3 * (m1 + m2) * m1 * m2 / C_SI**5
            * MYR_S / PARSEC_M**4)
    alpha = G_SI * rho * SOLAR_MASS_KG / PARSEC_M**3 / (velocity * 1000.0) * PARSEC_M * MYR_S
    polynomial = 1.0 + 73.0 / 24.0 * e2 + 37.0 / 96.0 * e2**2
    log_gw_rate = np.log(beta) - 4.0 * log_a - 3.5 * safe_log_j2 + np.log(polynomial)
    log_env_rate = np.log(alpha) + np.log(h) + log_a if alpha > 0.0 else np.full_like(a, -np.inf)
    log_rate = np.logaddexp(log_env_rate, log_gw_rate)
    # Combine exponentials before evaluating: no large intermediate 1/j².
    dy_du = (-2.0 * e * k * np.exp(log_env_rate - log_rate - safe_log_j2)
             + 19.0 / 6.0 * e2 * (1.0 + 121.0 / 304.0 * e2)
             / polynomial * np.exp(log_gw_rate - log_rate))
    answer = np.stack((np.exp(-log_rate), dy_du), axis=-1)
    return np.where(np.asarray(valid)[..., None], answer, np.nan)


def _evolve_adaptive_log_a(a, e, shell_state, shell_index, config):
    """逐轨道独立步长的向量化 Dormand-Prince 5(4)。

    使用安装版本 scipy.integrate.RK45 的 Butcher 系数和四次稠密输出系数。
    终点为 t=global_timestep 或 a=merger_radius_gm_c2*G*(m1+m2)/c²。
    后者是本项目的有限终止约定；不声称论文指定了 6GM/c²。
    Euler 的 local_timestep_myr 不控制此积分器，容差与 max_log_a_step 控制精度。
    """
    if config.primary_mass_msun != config.secondary_mass_msun:
        raise ValueError("当前 Sesana K 只实现 q=1。")
    rho = float(shell_state.environment_density_msun_pc3[shell_index])
    velocity = float(shell_state.velocity_dispersion_km_s[shell_index])
    ah = float(shell_state.hard_semimajor_axis_pc[shell_index])
    if not np.isfinite(rho) or rho < 0 or not np.isfinite(velocity) or velocity <= 0:
        raise ValueError("壳层密度和速度无效。")
    initial_a = a.copy()
    hard = a <= ah
    merged = np.zeros(a.size, dtype=bool)
    invalid = np.zeros(a.size, dtype=bool)
    k_above = hard & (e > 0.9)
    merger_step = np.full(a.size, -1, dtype=int)
    merger_time = np.full(a.size, np.nan)
    accepted_steps = np.zeros(a.size, dtype=int)
    cutoff = (config.merger_radius_gm_c2 * G_SI
              * (config.primary_mass_msun + config.secondary_mass_msun)
              * SOLAR_MASS_KG / C_SI**2 / PARSEC_M)
    u_stop = np.log(initial_a / cutoff)
    u = np.zeros(a.size)
    state = np.column_stack((np.zeros(a.size), np.log((1.0 - e) * (1.0 + e))))
    step_size = np.full(a.size, min(0.02, config.integration_max_log_a_step))
    pending = hard.copy()
    immediate = hard & (u_stop <= 0)
    merged[immediate] = True
    merger_time[immediate] = 0.0
    merger_step[immediate] = 0
    pending[immediate] = False
    duration = config.global_timestep_myr
    atol = np.array([config.integration_time_atol_myr, config.integration_log_j2_atol])

    for _ in range(config.integration_max_attempts):
        ix = np.flatnonzero(pending)
        if not ix.size:
            break
        uu, yy, aa = u[ix], state[ix], initial_a[ix]
        stages = np.empty((7, ix.size, 2))
        stages[0] = log_a_orbital_rhs(uu, yy, aa, rho, velocity, ah, config)
        # Limit the candidate near the time boundary; dense output locates its crossing.
        remaining = duration - yy[:, 0]
        hh = np.minimum(step_size[ix], u_stop[ix] - uu)
        hh = np.minimum(hh, 1.1 * remaining / np.maximum(stages[0, :, 0], 1e-300))
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            for stage in range(1, 6):
                trial = yy + hh[:, None] * np.einsum("s,snk->nk", RK45.A[stage, :stage], stages[:stage])
                stages[stage] = log_a_orbital_rhs(uu + RK45.C[stage] * hh, trial, aa, rho, velocity, ah, config)
            next_state = yy + hh[:, None] * np.einsum("s,snk->nk", RK45.B, stages[:6])
            stages[6] = log_a_orbital_rhs(uu + hh, next_state, aa, rho, velocity, ah, config)
            error = hh[:, None] * np.einsum("s,snk->nk", RK45.E, stages)
            scale = atol + config.integration_rtol * np.maximum(np.abs(yy), np.abs(next_state))
            error_norm = np.max(np.abs(error) / scale, axis=1)
        finite = np.all(np.isfinite(next_state), axis=1) & np.isfinite(error_norm)
        ok = finite & (error_norm <= 1) & (next_state[:, 1] <= 0) & (next_state[:, 0] >= yy[:, 0])
        factor = np.full(ix.size, 0.2)
        factor[finite] = np.clip(0.9 * np.maximum(error_norm[finite], 1e-16)**-0.2, 0.2, 5.0)
        factor[~ok] = np.minimum(factor[~ok], 0.5)
        step_size[ix] = np.minimum(config.integration_max_log_a_step, hh * factor)
        accepted = ix[ok]
        accepted_steps[accepted] += 1
        u[accepted] = uu[ok] + hh[ok]
        state[accepted] = next_state[ok]
        k_above[accepted] |= next_state[ok, 1] < np.log(0.19)

        at_time = ok & (next_state[:, 0] >= duration)
        if np.any(at_time):
            # Standard RK45 quartic dense output; bisection is vectorized per orbit.
            coeff = np.einsum("snk,sq->nkq", stages[:, at_time], RK45.P)
            lower = np.zeros(np.count_nonzero(at_time))
            upper = np.ones_like(lower)
            for _ in range(40):
                fraction = 0.5 * (lower + upper)
                powers = fraction[:, None] ** np.arange(1, 5)
                interpolated = yy[at_time] + hh[at_time, None] * np.einsum("nkq,nq->nk", coeff, powers)
                below = interpolated[:, 0] < duration
                lower = np.where(below, fraction, lower)
                upper = np.where(below, upper, fraction)
            ended = ix[at_time]
            u[ended] = uu[at_time] + hh[at_time] * fraction
            state[ended] = interpolated
            state[ended, 0] = duration
            pending[ended] = False

        at_merger = ok & ~at_time & (u[ix] >= u_stop[ix] - 1e-12)
        ended = ix[at_merger]
        merged[ended] = True
        merger_time[ended] = state[ended, 0]
        merger_step[ended] = accepted_steps[ended]
        pending[ended] = False
        stalled = pending[ix] & (step_size[ix] < 8 * np.spacing(np.maximum(1.0, uu)))
        invalid[ix[stalled]] = True
        pending[ix[stalled]] = False

    invalid[pending] = True
    a[hard] = initial_a[hard] * np.exp(-u[hard])
    a[merged] = 0.0
    e[hard] = np.sqrt(-np.expm1(state[hard, 1]))
    return BinaryEvolutionEvents(
        final_state=BinaryOrbitState(a, e, hard & ~merged & ~invalid),
        hard_mask=hard, merged_mask=merged, invalid_mask=invalid,
        k_above_calibration_mask=k_above, merger_step=merger_step,
        merger_time_myr=merger_time,
        completed_local_steps=int(accepted_steps.max()) if a.size else 0,
    )
