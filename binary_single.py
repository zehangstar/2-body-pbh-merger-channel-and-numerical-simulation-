"""Binary-single 的单时间片、单壳层批量轨道演化。

轨道方程消费 ``HaloShellState`` 中的环境密度、速度与硬双星边界；
``shell_density_msun_pc3`` 始终是未缩放的 NFW 局部密度。
这里仅返回原始 Monte Carlo 事件，不进行人口重加权或 HMF 积分。

论文原文印刷 Eq. (24) 对应 mu=2；Table III 诊断更接近 mu=0.5。
本项目以 Prada12-HMF、mu=0.5、ksi=1 为默认配置，任何表格匹配分支
必须在调用时显式指定，不能混称论文原始模型。
"""

from dataclasses import dataclass
import math

import numpy as np

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
    ``v_env = sqrt(mu * G * M(<r) / r)``。``ksi`` 与 PBH 的
    单体/双星人口比例是不同参数，不能互相代用。
    """

    halo_model: str = hc.DEFAULT_HALO_MODEL
    primary_mass_msun: float = 30.0
    secondary_mass_msun: float = 30.0
    pbh_fraction_total: float = 1.0
    fraction_single: float = 0.5
    fraction_binary: float = 0.5
    velocity_mu: float = 0.5
    ksi: float = 1.0
    global_timestep_myr: float = 200.0
    local_timestep_myr: float = 2.0
    start_redshift: float = 12.0
    sample_count_per_shell_and_step: int = 2_000_000

    def __post_init__(self):
        for name in (
            "primary_mass_msun",
            "secondary_mass_msun",
            "velocity_mu",
            "global_timestep_myr",
            "local_timestep_myr",
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
        steps = self.global_timestep_myr / self.local_timestep_myr
        if not math.isclose(steps, round(steps), rel_tol=0.0, abs_tol=1e-10):
            raise ValueError("全局时间步必须是局部时间步的整数倍。")

    @property
    def local_steps_per_global_step(self):
        """论文基准为 200/2=100 个 Euler 局部步。"""
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
    """一个冻结环境时间片的原始事件；没有物理双星数重加权。"""

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
    k = sesana_equal_mass_k(a, e, ah)
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
