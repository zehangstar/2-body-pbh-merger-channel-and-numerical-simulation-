"""Fig. 12 的未重加权 binary-single Monte Carlo 驱动层。

每个全局时间片独立从 Appendix D 联合分布抽取初态，逐壳、分块调用
``binary_single.evolve_binary_batch``，只累计原始样本并合数。论文第 IV 节
描述了每壳每步独立选取 ``N_sample``；没有说明未并合样本跨时间片的继承。
这里采用独立时间片口径，不声称它已由作者代码验证。

本模块不进行真实壳层双星数重加权、HMF 积分、平滑或绘图。
"""

from dataclasses import dataclass, replace
import math

import numpy as np
from scipy.optimize import brentq

import binary_single as bs
import haloconcentration as hc
import p_a_j


MPC_M = 1.0e6 * bs.PARSEC_M
MINIMUM_HALO_PBH_COUNT = 30.0


@dataclass(frozen=True)
class RawShellMergerHistory:
    """单条 ``M0`` 轨迹的 Fig. 12 原始计数，不含物理人口权重。

    ``time_edges_myr`` 是宇宙年龄；``*_per_step`` 的形状为
    ``(时间片数, 壳层数)``；
    ``cumulative_mergers`` 的形状为 ``(时间边界数, 壳层数)``，
    第一行是初始时刻的零并合数。
    """

    present_day_halo_mass_msun: float
    start_redshift_used: float
    configuration: bs.BinarySingleConfig
    halo_model: str
    sampler_name: str
    random_seed: int
    samples_per_shell_per_step: int
    chunk_size: int
    time_edges_myr: np.ndarray
    redshift_edges: np.ndarray
    local_timestep_per_step_myr: np.ndarray
    sampled_per_step: np.ndarray
    hard_per_step: np.ndarray
    mergers_per_step: np.ndarray
    invalid_per_step: np.ndarray
    k_above_calibration_per_step: np.ndarray
    cumulative_mergers: np.ndarray

    @property
    def shell_count(self):
        return self.mergers_per_step.shape[1]

    @property
    def time_step_count(self):
        return self.mergers_per_step.shape[0]

    @property
    def elapsed_time_edges_myr(self):
        """从模拟起点计时的 Fig. 12 横轴候选量。"""
        return self.time_edges_myr - self.time_edges_myr[0]


def cosmic_age_myr(redshift, cosmology):
    """平直物质+Lambda 背景下的宇宙年龄，单位 Myr。

    与 ``halostructure.rho_crit`` 采用相同的物质+Lambda 膨胀模型；
    在 z<=12 的本阶段时间网格中不加入辐射项。
    """
    z = np.asarray(redshift, dtype=float)
    if np.any(~np.isfinite(z)) or np.any(z < 0.0):
        raise ValueError("红移必须为有限非负数。")
    omega_m = float(cosmology["omega_m0"])
    omega_lambda = float(cosmology["omega_lambda0"])
    h = float(cosmology["h"])
    if omega_m <= 0.0 or omega_lambda <= 0.0 or h <= 0.0:
        raise ValueError("时间映射需要正的 Omega_m、Omega_Lambda 和 h。")
    if not math.isclose(omega_m + omega_lambda, 1.0, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError("当前解析时间映射只支持平直物质+Lambda 宇宙学。")
    hubble0_per_s = 100.0 * h * 1000.0 / MPC_M
    return (
        2.0
        / (3.0 * hubble0_per_s * np.sqrt(omega_lambda) * bs.MYR_S)
        * np.arcsinh(
            np.sqrt(omega_lambda / omega_m) / (1.0 + z) ** 1.5
        )
    )


def redshift_from_cosmic_age_myr(age_myr, cosmology):
    """``cosmic_age_myr`` 的解析逆映射。"""
    age = np.asarray(age_myr, dtype=float)
    if np.any(~np.isfinite(age)) or np.any(age <= 0.0):
        raise ValueError("宇宙年龄必须为有限正数。")
    omega_m = float(cosmology["omega_m0"])
    omega_lambda = float(cosmology["omega_lambda0"])
    h = float(cosmology["h"])
    if omega_m <= 0.0 or omega_lambda <= 0.0 or h <= 0.0:
        raise ValueError("时间映射需要正的 Omega_m、Omega_Lambda 和 h。")
    if not math.isclose(omega_m + omega_lambda, 1.0, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError("当前解析时间映射只支持平直物质+Lambda 宇宙学。")
    hubble0_per_s = 100.0 * h * 1000.0 / MPC_M
    argument = 1.5 * hubble0_per_s * np.sqrt(omega_lambda) * age * bs.MYR_S
    return (
        np.sqrt(omega_lambda / omega_m) / np.sinh(argument)
    ) ** (2.0 / 3.0) - 1.0


def global_time_grid(start_redshift, global_timestep_myr, cosmology):
    """从起始红移到 z=0 的时间边界；最后一片可以短于全局步。"""
    if not math.isfinite(global_timestep_myr) or global_timestep_myr <= 0.0:
        raise ValueError("全局时间步必须为有限正数。")
    start_age = float(cosmic_age_myr(start_redshift, cosmology))
    end_age = float(cosmic_age_myr(0.0, cosmology))
    if start_age >= end_age:
        raise ValueError("起始红移必须大于零。")
    full_step_count = int(np.floor((end_age - start_age) / global_timestep_myr))
    ages = start_age + global_timestep_myr * np.arange(full_step_count + 1)
    if end_age - ages[-1] > 1e-9:
        ages = np.append(ages, end_age)
    else:
        ages[-1] = end_age
    redshifts = redshift_from_cosmic_age_myr(ages, cosmology)
    redshifts[0] = start_redshift
    redshifts[-1] = 0.0
    return ages, redshifts


def start_redshift_with_minimum_pbh_count(present_day_halo_mass_msun, config):
    """将小 halo 的起点推迟到瞬时 PBH 数达到约 30 时。

    目标论文第 IV 节脚注：今天质量约不低于 ``5e4 M_sun`` 的 halo
    从 z=12 开始；更小 halo 从较低红移开始，使整个模拟期至少含约 30
    个 PBH。当前单色模型按 ``N_PBH(z)=f_PBH*M(z)/m`` 计数，其中
    ``m=primary_mass_msun``。这个计数也使低 ``f_PBH`` 的较大 halo
    在必要时推迟开始；它不涉及双星/单体比例或壳内人口重加权。
    """
    mass0 = float(present_day_halo_mass_msun)
    if not math.isfinite(mass0) or mass0 <= 0.0:
        raise ValueError("现今 halo 质量必须为有限正数。")
    if config.pbh_fraction_total <= 0.0:
        raise ValueError("要求至少 30 个 PBH 时，f_PBH 必须大于零。")
    minimum_halo_mass = (
        MINIMUM_HALO_PBH_COUNT
        * config.primary_mass_msun
        / config.pbh_fraction_total
    )
    if mass0 <= minimum_halo_mass:
        raise ValueError("该 halo 即使在 z=0 也没有超过约 30 个 PBH 的演化区间。")

    requested_start = config.start_redshift
    if requested_start <= 0.0:
        raise ValueError("模拟起始红移必须大于零。")

    def mass_above_threshold(z):
        return float(hc.mass_accretion_history(mass0, z, model=config.halo_model)) - minimum_halo_mass

    if mass_above_threshold(requested_start) >= 0.0:
        return requested_start
    return float(brentq(mass_above_threshold, 0.0, requested_start))


def run_raw_shell_monte_carlo(
    present_day_halo_mass_msun=1.15e12,
    config=None,
    orbital_distribution=None,
    random_seed=12345,
    chunk_size=100_000,
):
    """驱动 Fig. 12 的逐时间片、逐壳层、分块原始 Monte Carlo 计数。

    默认采样器是 ``p_a_j.AppendixDOrbitalDistribution`` 的严格联合分布，
    其 ``rho_eq`` 沿用该模块的显式默认值；可传入具有 ``sample(n,rng)``
    方法的分布替换。``random_seed`` 和 ``chunk_size`` 固定后可重复。

    每片使用片首红移更新一次 halo，内部环境冻结。最后不足一个全局步的
    时间片使用 ``ceil(剩余时间/局部步长)`` 个等长 Euler 步，以 z=0 收尾；
    这不是目标论文明确指定的终点算法。若请求的起点时 halo 内不足约
    30 个 PBH，则按 ``M(z)`` 把起点推迟；这里只保存计数，不保存个体轨道。
    """
    if config is None:
        config = bs.BinarySingleConfig()
    if not isinstance(random_seed, (int, np.integer)) or random_seed < 0:
        raise ValueError("random_seed 必须为非负整数。")
    if not isinstance(chunk_size, (int, np.integer)) or chunk_size <= 0:
        raise ValueError("chunk_size 必须为正整数。")
    n_sample = config.sample_count_per_shell_and_step
    if not isinstance(n_sample, (int, np.integer)) or n_sample <= 0:
        raise ValueError("sample_count_per_shell_and_step 必须为正整数。")
    if orbital_distribution is None:
        orbital_distribution = p_a_j.AppendixDOrbitalDistribution(
            pbh_mass_msun=config.primary_mass_msun,
            f_pbh=config.pbh_fraction_total,
        )
    if not callable(getattr(orbital_distribution, "sample", None)):
        raise TypeError("orbital_distribution 必须提供 sample(n, rng) 方法。")

    cosmology = hc.cosmology_for_model(config.halo_model)
    start_redshift = start_redshift_with_minimum_pbh_count(
        present_day_halo_mass_msun, config,
    )
    time_edges, redshift_edges = global_time_grid(
        start_redshift, config.global_timestep_myr, cosmology,
    )
    n_steps = time_edges.size - 1
    first_shell_state = bs.build_shell_environment(
        present_day_halo_mass_msun, redshift_edges[0], config,
    )
    n_shells = first_shell_state.shell_mass_msun.size
    shape = (n_steps, n_shells)
    sampled = np.zeros(shape, dtype=np.int64)
    hard = np.zeros(shape, dtype=np.int64)
    merged = np.zeros(shape, dtype=np.int64)
    invalid = np.zeros(shape, dtype=np.int64)
    k_above_calibration = np.zeros(shape, dtype=np.int64)
    actual_local_steps = np.zeros(n_steps, dtype=float)
    rng = np.random.default_rng(int(random_seed))

    for time_index in range(n_steps):
        duration = float(time_edges[time_index + 1] - time_edges[time_index])
        if math.isclose(duration, config.global_timestep_myr, rel_tol=0.0, abs_tol=1e-9):
            step_config = config
        else:
            local_count = max(1, math.ceil(duration / config.local_timestep_myr))
            step_config = replace(
                config,
                global_timestep_myr=duration,
                local_timestep_myr=duration / local_count,
            )
        actual_local_steps[time_index] = step_config.local_timestep_myr
        shell_state = (
            first_shell_state if time_index == 0
            else bs.build_shell_environment(
                present_day_halo_mass_msun, redshift_edges[time_index], config,
            )
        )
        if shell_state.shell_mass_msun.size != n_shells:
            raise ValueError("同一 M0 轨迹的壳层数不能随红移改变。")

        for shell_index in range(n_shells):
            remaining = n_sample
            while remaining:
                draw_count = min(remaining, chunk_size)
                orbits = orbital_distribution.sample(draw_count, rng)
                events = bs.evolve_binary_batch(
                    orbits.semi_major_axis_pc,
                    orbits.eccentricity,
                    shell_state,
                    shell_index,
                    step_config,
                )
                if events.sample_count != draw_count:
                    raise ValueError("采样器返回数量与请求的分块样本数不一致。")
                sampled[time_index, shell_index] += events.sample_count
                hard[time_index, shell_index] += events.hard_count
                merged[time_index, shell_index] += events.merger_count
                invalid[time_index, shell_index] += events.invalid_count
                k_above_calibration[time_index, shell_index] += int(
                    np.count_nonzero(events.k_above_calibration_mask)
                )
                remaining -= draw_count

    cumulative = np.vstack((
        np.zeros((1, n_shells), dtype=np.int64),
        np.cumsum(merged, axis=0),
    ))
    return RawShellMergerHistory(
        present_day_halo_mass_msun=float(present_day_halo_mass_msun),
        start_redshift_used=start_redshift,
        configuration=config,
        halo_model=config.halo_model,
        sampler_name=type(orbital_distribution).__name__,
        random_seed=int(random_seed),
        samples_per_shell_per_step=int(n_sample),
        chunk_size=int(chunk_size),
        time_edges_myr=time_edges,
        redshift_edges=redshift_edges,
        local_timestep_per_step_myr=actual_local_steps,
        sampled_per_step=sampled,
        hard_per_step=hard,
        mergers_per_step=merged,
        invalid_per_step=invalid,
        k_above_calibration_per_step=k_above_calibration,
        cumulative_mergers=cumulative,
    )
