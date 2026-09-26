"""Binary-single Monte Carlo：方案 A 独立重抽和方案 B cohort 继承。

每个全局时间片独立从 Appendix D 联合分布抽取初态，逐壳、分块调用
``binary_single.evolve_binary_batch``，只累计原始样本并合数。论文第 IV 节
描述了每壳每步独立选取 ``N_sample``；没有说明未并合样本跨时间片的继承。
这里采用独立时间片口径，不声称它已由作者代码验证。

run_raw_shell_monte_carlo 保留方案 A 的原始样本计数。
run_cohort_shell_monte_carlo 按新增壳层质量赋予权重，保存和推进存活轨道。
两者都不做 HMF 积分、平滑或绘图。
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
    local_timestep_per_step_myr 只记录 Euler 实际步长；自适应分支为 NaN，
    精度设置保存在 configuration，不能把遗留的 2 Myr 参数当作其实际步长。
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
    progress_callback=None,
):
    """驱动 Fig. 12 的逐时间片、逐壳层、分块原始 Monte Carlo 计数。

    默认采样器是 ``p_a_j.AppendixDOrbitalDistribution`` 的严格联合分布，
    其 ``rho_eq`` 沿用该模块的显式默认值；可传入具有 ``sample(n,rng)``
    方法的分布替换。``random_seed`` 和 ``chunk_size`` 固定后可重复。
    可选 progress_callback(completed_steps, total_steps) 在每片完成后调用。

    每片使用片首红移更新一次 halo，内部环境冻结。最后不足一个全局步的
    时间片在 Euler 分支使用 ``ceil(剩余时间/局部步长)`` 个等长步；
    自适应分支直接以实际剩余时长终止，两者均以 z=0 收尾；
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
        elif config.integration_method == "euler":
            local_count = max(1, math.ceil(duration / config.local_timestep_myr))
            step_config = replace(
                config,
                global_timestep_myr=duration,
                local_timestep_myr=duration / local_count,
            )
        else:
            step_config = replace(config, global_timestep_myr=duration)
        actual_local_steps[time_index] = (
            step_config.local_timestep_myr
            if config.integration_method == "euler" else np.nan
        )
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
        if progress_callback is not None:
            progress_callback(time_index + 1, n_steps)

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


@dataclass(frozen=True)
class CohortPopulationConfig:
    """cohort 人口模型的显式选择，不代表原论文未公开的实现。

    entry_orbit_model='gw_aged' 将原初轨道在晕外以 GW 演化至入晕时刻；
    'pristine' 仅用于隔离人口继承影响的对照。formation_age_myr=0 是
    等时形成于宇宙早期的时间零近似，可替换为指定形成宇宙年龄。
    新增人口按每壳质量增量注入，个体保持壳编号；因此只支持每壳质量
    单调增长的轨迹。不对负增量取绝对值或静默补回人口。
    """
    initial_samples_per_shell: int = 5000
    accretion_samples_per_shell: int = 256
    entry_orbit_model: str = 'gw_aged'
    formation_age_myr: float = 0.0

    def __post_init__(self):
        for n in (self.initial_samples_per_shell, self.accretion_samples_per_shell):
            if not isinstance(n, (int, np.integer)) or n <= 0:
                raise ValueError('cohort 样本数必须是正整数。')
        if self.entry_orbit_model not in ('gw_aged', 'pristine'):
            raise ValueError('entry_orbit_model 必须是 gw_aged 或 pristine。')
        if not math.isfinite(self.formation_age_myr) or self.formation_age_myr < 0:
            raise ValueError('formation_age_myr 必须是有限非负数。')


@dataclass(frozen=True)
class CohortOrbitalPopulation:
    """单壳存活者；entry_boundary 保留入晕 cohort 身份。"""
    semi_major_axis_pc: np.ndarray
    eccentricity: np.ndarray
    weight: np.ndarray
    entry_boundary: np.ndarray


@dataclass(frozen=True)
class CohortShellMergerHistory:
    present_day_halo_mass_msun: float
    configuration: bs.BinarySingleConfig
    population_configuration: CohortPopulationConfig
    random_seed: int
    time_edges_myr: np.ndarray
    redshift_edges: np.ndarray
    shell_mass_at_boundary_msun: np.ndarray
    drawn_per_boundary: np.ndarray
    supplied_binary_weight_per_boundary: np.ndarray
    outside_merged_weight_per_boundary: np.ndarray
    entered_alive_weight_per_boundary: np.ndarray
    surviving_weight_at_boundary: np.ndarray
    raw_mergers_per_step: np.ndarray
    merged_weight_per_step: np.ndarray
    soft_gw_merged_weight_per_step: np.ndarray
    final_populations: tuple

    @property
    def elapsed_time_edges_myr(self):
        return self.time_edges_myr - self.time_edges_myr[0]

    @property
    def rate_shell_per_year(self):
        return self.merged_weight_per_step / (np.diff(self.time_edges_myr)[:, None] * 1e6)

    @property
    def cumulative_merger_weight(self):
        return np.vstack((np.zeros(self.merged_weight_per_step.shape[1]),
                          np.cumsum(self.merged_weight_per_step, axis=0)))

    @property
    def population_balance_residual(self):
        """逐边界：累计供给 = 晕外已并合 + 晕内已并合 + 当前存活。"""
        return (np.cumsum(self.supplied_binary_weight_per_boundary, axis=0)
                - np.cumsum(self.outside_merged_weight_per_boundary, axis=0)
                - self.cumulative_merger_weight - self.surviving_weight_at_boundary)


def run_cohort_shell_monte_carlo(
    present_day_halo_mass_msun=1.15e12,
    config=None,
    population_config=None,
    orbital_distribution=None,
    random_seed=12345,
    progress_callback=None,
    checkpoint_path=None,
    resume=False,
):
    """方案 B：带入晕 cohort 和质量权重的人口守恒轨道演化。

    每个边界 k 将 f_PBH*f_binary*Delta M_shell/(m1+m2) 的原初双星
    系统权重注入该壳。每壳新注入 n 个独立样本，每个权重 Delta N/n。
    gw_aged 入口会先扣除晕外已经并合的权重，并保留入晕前 GW 演化状态。
    晕内硬双星用完整轨道方程；软双星仍以 GW 演化。每片永久移除并合者，
    存活者保留 a,e,weight,entry_boundary；绝不乘回完整壳层人口。

    环境冻结于片首，新增质量在时间边界加入；最终 z=0 边界也登记新质量，
    但不给它额外演化时长。final_populations 包含这些刚进入的终点存活者。
    保持壳编号、按各壳质量增量供给、早期共同形成是显式模型闭合条件。
    不含反冲、软双星离解、再形成、双星合并后代再入或壳间迁移。

    可选 checkpoint_path 在每个完整时间边界后保存人口、账本和 RNG 状态；
    resume=True 从该文件继续（文件尚不存在则从头开始）。配置、源码或
    严格采样器参数改变时拒绝续算；不改变无检查点调用的抽样顺序。
    """
    config = bs.BinarySingleConfig() if config is None else config
    population_config = CohortPopulationConfig() if population_config is None else population_config
    if config.integration_method != 'adaptive_log_a':
        raise ValueError('cohort 的任意入晕年龄目前需要 adaptive_log_a。')
    if not isinstance(random_seed, (int, np.integer)) or random_seed < 0:
        raise ValueError('random_seed 必须是非负整数。')
    distribution = orbital_distribution
    if distribution is None:
        distribution = p_a_j.AppendixDOrbitalDistribution(
            pbh_mass_msun=config.primary_mass_msun, f_pbh=config.pbh_fraction_total)
    if not callable(getattr(distribution, 'sample', None)):
        raise TypeError('orbital_distribution 必须提供 sample(n,rng)。')
    cosmology = hc.cosmology_for_model(config.halo_model)
    z0 = start_redshift_with_minimum_pbh_count(present_day_halo_mass_msun, config)
    times, redshifts = global_time_grid(z0, config.global_timestep_myr, cosmology)
    if population_config.formation_age_myr > times[0]:
        raise ValueError('形成时间不能晚于首次入晕时间。')
    shells = [bs.build_shell_environment(present_day_halo_mass_msun, float(z), config) for z in redshifts]
    masses = np.array([s.shell_mass_msun for s in shells])
    increments = np.diff(np.vstack((np.zeros(masses.shape[1]), masses)), axis=0)
    tolerance = 1e-12 * max(float(np.max(masses)), 1.0)
    if np.any(increments < -tolerance):
        raise ValueError('壳质量有负增量：需要显式壳间输运模型，不能使用当前 cohort 入口。')
    increments = np.maximum(increments, 0.0)  # 仅容忍上方检查通过的舍入误差。
    supplied = (config.pbh_fraction_total * config.fraction_binary * increments
                / (config.primary_mass_msun + config.secondary_mass_msun))
    nsteps, nshell = len(times)-1, masses.shape[1]
    drawn = np.zeros_like(masses, dtype=np.int64)
    outside = np.zeros_like(masses)
    entered = np.zeros_like(masses)
    surviving = np.zeros_like(masses)
    raw_mergers = np.zeros((nsteps, nshell), dtype=np.int64)
    merged_weight = np.zeros((nsteps, nshell))
    soft_merged_weight = np.zeros_like(merged_weight)
    populations = [CohortOrbitalPopulation(np.empty(0),np.empty(0),np.empty(0),np.empty(0,dtype=int)) for _ in range(nshell)]
    rng = np.random.default_rng(int(random_seed))
    first_boundary = 0
    checkpoint = None
    if resume and checkpoint_path is None:
        raise ValueError('resume=True 需要 checkpoint_path。')
    if checkpoint_path is not None:
        import hashlib
        import json
        from dataclasses import asdict
        from pathlib import Path
        if type(distribution) is not p_a_j.AppendixDOrbitalDistribution:
            raise ValueError('检查点目前只支持显式参数可记录的严格 Appendix D 采样器。')
        checkpoint = Path(checkpoint_path)
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        identity = dict(schema=1, mass_msun=float(present_day_halo_mass_msun),
            config=asdict(config), population=asdict(population_config),
            distribution=asdict(distribution), seed=int(random_seed),
            sources={Path(p).name: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                     for p in (__file__, bs.__file__, hc.__file__, bs.hs.__file__, p_a_j.__file__)})
        ledger_names = ('drawn', 'outside', 'entered', 'surviving', 'raw_mergers',
                        'merged_weight', 'soft_merged_weight')
        ledgers = (drawn, outside, entered, surviving, raw_mergers,
                   merged_weight, soft_merged_weight)
        if checkpoint.exists():
            if not resume:
                raise FileExistsError('检查点已存在；请使用 resume=True 或新的输出目录。')
            with np.load(checkpoint, allow_pickle=False) as saved:
                if json.loads(str(saved['identity'])) != identity:
                    raise ValueError('检查点配置/采样器/源码与当前运行不一致，不能续算。')
                if not np.allclose(saved['time_edges_myr'], times, rtol=0, atol=1e-9):
                    raise ValueError('检查点时间网格不一致。')
                first_boundary = int(saved['next_boundary'])
                if not 0 <= first_boundary <= nsteps + 1:
                    raise ValueError('检查点边界索引无效。')
                for name, target in zip(ledger_names, ledgers):
                    target[:] = saved[name]
                rng.bit_generator.state = json.loads(str(saved['rng_state']))
                populations = [CohortOrbitalPopulation(*[
                    saved[f'population_{i}_{field}'].copy()
                    for field in ('a', 'e', 'weight', 'entry')]) for i in range(nshell)]

    for k in range(first_boundary, nsteps+1):
        shell = shells[k]
        for i in range(nshell):
            if supplied[k,i] > 0:
                n = (population_config.initial_samples_per_shell if k == 0
                     else population_config.accretion_samples_per_shell)
                sample = distribution.sample(n, rng)
                a, e = sample.semi_major_axis_pc.copy(), sample.eccentricity.copy()
                if a.size != n:
                    raise ValueError('cohort 入口采样数量不一致。')
                weight = supplied[k,i]/n
                drawn[k,i] = n
                age = times[k] - population_config.formation_age_myr
                if population_config.entry_orbit_model == 'gw_aged' and age > 0:
                    isolated = replace(shell,
                        environment_density_msun_pc3=np.zeros(nshell),
                        hard_semimajor_axis_pc=np.full(nshell, np.max(a)))
                    aged = bs.evolve_binary_batch(a,e,isolated,i,replace(config,
                        global_timestep_myr=float(age),eccentricity_growth_model='zero'))
                    if aged.invalid_count:
                        raise RuntimeError(f'晕外演化失效：边界 {k}，壳 {i+1}，{aged.invalid_count} 条。')
                    outside[k,i] = weight*aged.merger_count
                    keep = ~aged.merged_mask
                    a, e = aged.final_state.semi_major_axis_pc[keep], aged.final_state.eccentricity[keep]
                entered[k,i] = weight*a.size
                old = populations[i]
                populations[i] = CohortOrbitalPopulation(
                    np.r_[old.semi_major_axis_pc,a], np.r_[old.eccentricity,e],
                    np.r_[old.weight,np.full(a.size,weight)],
                    np.r_[old.entry_boundary,np.full(a.size,k,dtype=int)])
            population = populations[i]
            surviving[k,i] = population.weight.sum()
            if k == nsteps or not population.weight.size:
                continue
            step_config = replace(config,global_timestep_myr=float(times[k+1]-times[k]))
            event = bs.evolve_binary_population_batch(population.semi_major_axis_pc,
                population.eccentricity,shell,i,step_config)
            if event.invalid_count:
                raise RuntimeError(f'晕内演化失效：时间片 {k}，壳 {i+1}，{event.invalid_count} 条。')
            raw_mergers[k,i] = event.merger_count
            merged_weight[k,i] = population.weight[event.merged_mask].sum()
            soft_merged_weight[k,i] = population.weight[event.merged_mask & ~event.hard_mask].sum()
            keep = ~event.merged_mask
            populations[i] = CohortOrbitalPopulation(
                event.final_state.semi_major_axis_pc[keep],event.final_state.eccentricity[keep],
                population.weight[keep],population.entry_boundary[keep])
        if checkpoint is not None:
            payload = dict(zip(ledger_names, ledgers))
            payload.update(identity=np.array(json.dumps(identity, sort_keys=True)),
                           rng_state=np.array(json.dumps(rng.bit_generator.state)),
                           next_boundary=np.array(k+1), time_edges_myr=times)
            for i, population in enumerate(populations):
                for field, value in zip(('a', 'e', 'weight', 'entry'),
                        (population.semi_major_axis_pc, population.eccentricity,
                         population.weight, population.entry_boundary)):
                    payload[f'population_{i}_{field}'] = value
            temporary = checkpoint.with_suffix(checkpoint.suffix + '.partial')
            with temporary.open('wb') as stream:
                np.savez_compressed(stream, **payload)
            temporary.replace(checkpoint)
        if progress_callback is not None and k < nsteps:
            progress_callback(k+1,nsteps)
    return CohortShellMergerHistory(
        present_day_halo_mass_msun=float(present_day_halo_mass_msun),
        configuration=config,population_configuration=population_config,
        random_seed=int(random_seed),time_edges_myr=times,redshift_edges=redshifts,
        shell_mass_at_boundary_msun=masses,drawn_per_boundary=drawn,
        supplied_binary_weight_per_boundary=supplied,outside_merged_weight_per_boundary=outside,
        entered_alive_weight_per_boundary=entered,surviving_weight_at_boundary=surviving,
        raw_mergers_per_step=raw_mergers,merged_weight_per_step=merged_weight,
        soft_gw_merged_weight_per_step=soft_merged_weight,final_populations=tuple(populations))
