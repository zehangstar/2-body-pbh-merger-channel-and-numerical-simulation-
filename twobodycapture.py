"""暗物质晕中的 PBH 两体引力波捕获。

本模块实现目标论文 Eq. (11) 的捕获截面、速度平均 ``<Sigma v>``、
Appendix A14 对任意 PBH 质量分布的每晕捕获率，以及 Eq. (18) 对晕质量函数
积分得到的总共动捕获率。单色和连续质量谱使用同一个每晕率函数，区别只在
传入的质量积分节点和权重。

速度分布沿用 ``halostructure.py`` 中归一化的三维各向同性密度 ``f_3(v)``：

    4*pi*integral v^2*f_3(v) dv = 1.

因此速度平均显式使用球壳测度 ``4*pi*v^2 dv``，不会重复或遗漏 Jacobian。
这里没有自动自检、断言或独立测试模块。
"""

import numpy as np
from astropy.cosmology import Planck18
from hmf import MassFunction
from scipy.integrate import quad

import haloconcentration as hc
import halostructure as hs


SPEED_OF_LIGHT_KM_S = 299792.458
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
KILOMETERS_PER_KPC = 3.0856775814913673e16
KM_S_TO_KPC_PER_YEAR = SECONDS_PER_YEAR / KILOMETERS_PER_KPC

# 论文只说明晕质量函数由 HMFcalc 计算，并未列出 HMFcalc 的完整配置。
# HMFcalc 是网页工具；这里安装并导入的 ``hmf`` 是它实际使用的 Python 后端。
# 为避免结果随包的默认值变化，先把本阶段采用的 Planck18 功率谱参数明确写出。
# 这些数值是当前 hmf 3.6.3 的 Planck18 默认值，后续若找到作者配置应在此替换。
HMF_SIGMA_8 = 0.8159
HMF_SPECTRAL_INDEX = 0.9667
MPC3_PER_GPC3 = 1.0e9

def gw_capture_cross_section_kpc2(
    mass1_msun,
    mass2_msun,
    relative_speed_km_s,
):
    """计算论文 Eq. (11) 的两体引力波捕获截面，单位为 kpc^2。

    输入的两个 PBH 质量使用 M_sun，相对速度使用 km/s。这里的 ``sigma``
    是相互作用截面，不是密度涨落 ``sigma(M,z)``，也不是速度弥散。
    """
    mass1_msun = np.asarray(mass1_msun, dtype=float)
    mass2_msun = np.asarray(mass2_msun, dtype=float)
    relative_speed_km_s = np.asarray(relative_speed_km_s, dtype=float)

    if np.any(mass1_msun <= 0.0) or np.any(mass2_msun <= 0.0):
        raise ValueError("两个 PBH 质量都必须为正。")
    if np.any(relative_speed_km_s <= 0.0):
        raise ValueError("相对速度必须为正。")

    numerical_factor = 2.0 * np.pi * (
        85.0 * np.pi / (6.0 * np.sqrt(2.0))
    ) ** (2.0 / 7.0)
    mass_factor = (
        (mass1_msun + mass2_msun) ** (10.0 / 7.0)
        * mass1_msun ** (2.0 / 7.0)
        * mass2_msun ** (2.0 / 7.0)
    )
    return (
        numerical_factor
        * hs.GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN**2
        * mass_factor
        / (
            SPEED_OF_LIGHT_KM_S ** (10.0 / 7.0)
            * relative_speed_km_s ** (18.0 / 7.0)
        )
    )


def velocity_averaged_capture_coefficient_kpc3_per_year(
    mass1_msun,
    mass2_msun,
    velocity_dispersion_km_s,
    cutoff_velocity_km_s,
):
    """计算 ``<Sigma_capture*v>``，返回单位 kpc^3/yr。

    直接使用无歧义的定义

        <Sigma v> = 4*pi*integral v^2*f_3(v)*Sigma(v)*v dv.

    将 Eq. (11) 的 ``v^(-18/7)`` 与测度中的 ``v^3`` 合并后，被积函数
    在低速端正比于 ``v^(3/7)``，因此不会在 v=0 发散。
    """
    velocity_dispersion_km_s = float(velocity_dispersion_km_s)
    cutoff_velocity_km_s = float(cutoff_velocity_km_s)
    if velocity_dispersion_km_s <= 0.0 or cutoff_velocity_km_s <= 0.0:
        raise ValueError("速度弥散和截止速度必须为正。")

    def velocity_integrand(relative_speed_km_s):
        velocity_density = float(
            hs.truncated_maxwell_velocity_density(
                relative_speed_km_s,
                velocity_dispersion_km_s,
                cutoff_velocity_km_s,
            )
        )
        return (
            4.0
            * np.pi
            * relative_speed_km_s ** (3.0 / 7.0)
            * velocity_density
        )

    velocity_moment, _ = quad(
        velocity_integrand,
        0.0,
        cutoff_velocity_km_s,
        epsabs=0.0,
        epsrel=1.0e-9,
        limit=200,
    )

    mass1_msun = np.asarray(mass1_msun, dtype=float)
    mass2_msun = np.asarray(mass2_msun, dtype=float)
    if np.any(mass1_msun <= 0.0) or np.any(mass2_msun <= 0.0):
        raise ValueError("两个 PBH 质量都必须为正。")

    numerical_factor = 2.0 * np.pi * (
        85.0 * np.pi / (6.0 * np.sqrt(2.0))
    ) ** (2.0 / 7.0)
    mass_factor = (
        (mass1_msun + mass2_msun) ** (10.0 / 7.0)
        * mass1_msun ** (2.0 / 7.0)
        * mass2_msun ** (2.0 / 7.0)
    )

    # 上式先得到 kpc^2*(km/s)，再把 km/s 换成 kpc/yr。
    return (
        numerical_factor
        * hs.GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN**2
        * mass_factor
        * velocity_moment
        / SPEED_OF_LIGHT_KM_S ** (10.0 / 7.0)
        * KM_S_TO_KPC_PER_YEAR
    )


def capture_rate_per_halo_a14_per_year(
    halo_mass_msun,
    z,
    mass_nodes_msun,
    mass_probability_weights,
    concentration_model="ludlow16",
    f_pbh=1.0,
    concentration_override=None,
):
    """用 Appendix A14 计算任意 PBH 质量分布的每晕捕获率，单位 yr^-1。

    ``mass_nodes_msun`` 与 ``mass_probability_weights`` 共同表示质量测度：

        integral psi(m) F(m) dm ~= sum_i weight_i F(m_i).

    单色谱传入一个节点和权重 1；连续谱传入数值积分节点与 ``psi(m)dm``
    权重。两者调用完全相同的本函数，不存在单色专用每晕率公式。

    ``f_pbh`` 是参与直接捕获的有效 PBH 质量占暗物质的比例。论文基准在
    Eq. (12) 附近取总 PBH 丰度为 1；若只允许部分 PBH 参与，应显式传入相应值，
    结果按 ``f_pbh**2`` 缩放。

    ``concentration_override`` 默认是 ``None``，此时函数按模型名自行计算
    ``C(M,z)``。Fig. 8 的 Prada12 复现需要用 HMF 功率谱计算真实的
    ``sigma(M,z)``，调用者可把已算好的浓度传入这里；NFW、速度分布和捕获率
    的其余计算仍走同一条路径。
    """
    halo_mass_msun = float(halo_mass_msun)
    z = float(z)
    mass_nodes_msun = np.asarray(mass_nodes_msun, dtype=float)
    mass_probability_weights = np.asarray(
        mass_probability_weights, dtype=float
    )

    if halo_mass_msun <= 0.0:
        raise ValueError("晕质量必须为正。")
    if mass_nodes_msun.ndim != 1 or mass_probability_weights.ndim != 1:
        raise ValueError("质量节点和权重必须是一维数组。")
    if mass_nodes_msun.shape != mass_probability_weights.shape:
        raise ValueError("质量节点和权重必须具有相同长度。")
    if np.any(mass_nodes_msun <= 0.0) or np.any(
        mass_probability_weights < 0.0
    ):
        raise ValueError("质量节点必须为正，质量权重不能为负。")
    if f_pbh < 0.0:
        raise ValueError("f_pbh 不能为负。")

    model_name = concentration_model.lower()
    if model_name == "ludlow16":
        cosmology = hc.LUDLOW_COSMOLOGY
    elif model_name in {"prada12", "prada12_hmf_sigma"}:
        cosmology = hc.PRADA_COSMOLOGY
    else:
        raise ValueError(
            "concentration_model 必须是 'ludlow16'、'prada12' 或 "
            "'prada12_hmf_sigma'。"
        )

    if concentration_override is None:
        concentration = float(
            hc.concentration_model(
                halo_mass_msun,
                z,
                model_name,
                cap_prada=(
                    model_name in {"prada12", "prada12_hmf_sigma"}
                ),
            )
        )
    else:
        concentration = float(concentration_override)
        if not np.isfinite(concentration) or concentration <= 0.0:
            raise ValueError("外部传入的 concentration 必须是有限正数。")
    critical_density_msun_kpc3 = float(hs.rho_crit(z, cosmology))
    virial_radius_kpc = (
        3.0
        * halo_mass_msun
        / (4.0 * np.pi * 200.0 * critical_density_msun_kpc3)
    ) ** (1.0 / 3.0)
    scale_radius_kpc = virial_radius_kpc / concentration
    scale_density_msun_kpc3 = float(
        hs.rho_s(concentration, z, cosmology=cosmology)
    )
    velocity_dispersion_km_s = float(
        hs.nfw_characteristic_velocity(
            scale_radius_kpc,
            scale_density_msun_kpc3,
        )
    )
    cutoff_velocity_km_s = float(
        hs.virial_cutoff_velocity(halo_mass_msun, virial_radius_kpc)
    )

    mass1_msun = mass_nodes_msun[:, np.newaxis]
    mass2_msun = mass_nodes_msun[np.newaxis, :]
    pair_probability_weights = (
        mass_probability_weights[:, np.newaxis]
        * mass_probability_weights[np.newaxis, :]
    )
    sigma_v_kpc3_per_year = (
        velocity_averaged_capture_coefficient_kpc3_per_year(
            mass1_msun,
            mass2_msun,
            velocity_dispersion_km_s,
            cutoff_velocity_km_s,
        )
    )
    mass_averaged_pair_coefficient = np.sum(
        pair_probability_weights
        * sigma_v_kpc3_per_year
        / (mass1_msun * mass2_msun)
    )

    # 这是 4*pi*integral r^2*rho_NFW(r)^2 dr 的解析结果。
    concentration_factor = 1.0 - 1.0 / (1.0 + concentration) ** 3
    density_squared_volume_integral = (
        halo_mass_msun**2
        * concentration_factor
        / (
            12.0
            * np.pi
            * scale_radius_kpc**3
            * hc.nfw_g(concentration) ** 2
        )
    )

    # 完整 (m1,m2) 正方形会把无序质量对计算两次，所以保留 Appendix A 的 1/2。
    return float(
        0.5
        * f_pbh**2
        * density_squared_volume_integral
        * mass_averaged_pair_coefficient
    )


def press_schechter_halo_mass_function_hmfcalc(
    z,
    minimum_halo_mass_msun=1.0e3,
    maximum_halo_mass_msun=1.0e15,
    dlog10m=0.1,
):
    """用 HMFcalc 的 ``hmf`` 后端计算 Press-Schechter 晕质量函数。

    论文 Eq. (18) 需要 ``dn/dM``，但跨越许多个数量级时，更适合使用与它
    完全等价的对数形式

        (dn/dM) dM = (dn/dlnM) dlnM.

    因此本函数返回两个一维数组：

    ``halo_masses_msun``
        晕的物理质量 M，单位 M_sun。

    ``dndlnm_mpc3``
        每单位 ``ln(M)`` 的共动晕数密度，单位 Mpc^-3。

    ``hmf`` 内部的质量坐标是 ``M_sun/h``，``dndlnm`` 的原始单位是
    ``h^3 Mpc^-3``。下面分别除以 ``h`` 和乘以 ``h^3``，把公开接口统一为
    物理 ``M_sun`` 与 ``Mpc^-3``。自然对数在乘以常数 ``h`` 前后微分相同，
    即 ``dln(M/h)=dln(M)``，所以不再产生额外 Jacobian。

    注意：论文没有给出作者使用的完整 HMFcalc 配置。本阶段明确采用
    Planck18、CAMB 转移函数、Press-Schechter 拟合、``delta_c=1.686``，以及
    模块顶部列出的 ``sigma_8`` 和谱指数。这里是可追踪的复现假设，不应被
    误认为论文已经逐项确认的参数。
    """
    z = float(z)
    minimum_halo_mass_msun = float(minimum_halo_mass_msun)
    maximum_halo_mass_msun = float(maximum_halo_mass_msun)
    dlog10m = float(dlog10m)

    if z < 0.0:
        raise ValueError("红移 z 不能为负。")
    if minimum_halo_mass_msun <= 0.0:
        raise ValueError("最小晕质量必须为正。")
    if maximum_halo_mass_msun <= minimum_halo_mass_msun:
        raise ValueError("最大晕质量必须大于最小晕质量。")
    if dlog10m <= 0.0:
        raise ValueError("dlog10m 必须为正。")

    hubble_h = float(Planck18.h)

    # MassFunction 的 Mmin/Mmax 是 log10(M/[M_sun/h])。
    # Mmax 再增加半个步长，使目标上边界能够进入 hmf 的左闭右开质量网格。
    mass_function = MassFunction(
        Mmin=np.log10(minimum_halo_mass_msun * hubble_h),
        Mmax=(
            np.log10(maximum_halo_mass_msun * hubble_h)
            + 0.5 * dlog10m
        ),
        dlog10m=dlog10m,
        z=z,
        hmf_model="PS",
        delta_c=1.686,
        cosmo_model=Planck18,
        sigma_8=HMF_SIGMA_8,
        n=HMF_SPECTRAL_INDEX,
        transfer_model="CAMB",
        transfer_params={"extrapolate_with_eh": True},
    )

    halo_masses_msun = np.asarray(mass_function.m, dtype=float) / hubble_h
    dndlnm_mpc3 = (
        np.asarray(mass_function.dndlnm, dtype=float) * hubble_h**3
    )

    # 浮点运算有时会把理论上的 10^3 写成 999.9999999999999，因此边界筛选
    # 留出极小的相对容差；这不是物理参数调整。
    boundary_tolerance = 1.0e-12
    inside_requested_range = (
        halo_masses_msun
        >= minimum_halo_mass_msun * (1.0 - boundary_tolerance)
    ) & (
        halo_masses_msun
        <= maximum_halo_mass_msun * (1.0 + boundary_tolerance)
    )

    return (
        halo_masses_msun[inside_requested_range],
        dndlnm_mpc3[inside_requested_range],
    )


def comoving_capture_rate_with_history_gpc3_per_year(
    z,
    mass_nodes_msun,
    mass_probability_weights,
    mass_history_function,
    concentration_model="ludlow16",
    f_pbh=1.0,
    minimum_present_halo_mass_msun=1.0e3,
    maximum_present_halo_mass_msun=1.0e15,
    halo_track_point_count=50,
    hmf_dlog10m=0.05,
):
    """沿给定质量吸积史计算总共动捕获率。

    ``mass_history_function`` 必须接受
    ``(present_day_masses, z, concentration_model)``，并返回每个今天质量
    ``M0`` 在红移 ``z`` 时对应的实际晕质量 ``M(z;M0)``。Fig. 8 专属的
    吸积史数据和插值函数保存在绘图 Notebook，而不是本物理模块中。

    原始矢量图的曲线可由以下过程重现到约 10%：

    1. 在今天的质量 ``M0=10^3--10^15 M_sun`` 上取 50 个对数等距晕；
    2. 沿各自质量吸积史得到非均匀的 ``M(z;M0)``，并计算每晕率；
    3. 在 ``min[M(z)]--max[M(z)]`` 之间另建均匀对数 HMF 网格；
    4. 按数组下标把每晕率与 HMF 值配对，再对 ``ln(M0)`` 积分。

    第 3--4 步是从论文原始矢量数据反推的“作图数值约定”，论文正文没有公开
    对应代码。它不同于严格地在每个实际 ``M(z;M0)`` 上评价 HMF，因此不能把
    二者的差异解释成新的物理效应。保留现有
    :func:`comoving_capture_rate_eq18_gpc3_per_year`，就是为了让严格的当前质量
    积分与论文图复现不被混成同一个接口。
    """
    z = float(z)
    minimum_present_halo_mass_msun = float(minimum_present_halo_mass_msun)
    maximum_present_halo_mass_msun = float(maximum_present_halo_mass_msun)
    halo_track_point_count = int(halo_track_point_count)

    if z < 0.0:
        raise ValueError("红移 z 不能为负。")
    if minimum_present_halo_mass_msun <= 0.0:
        raise ValueError("今天的最小晕质量必须为正。")
    if maximum_present_halo_mass_msun <= minimum_present_halo_mass_msun:
        raise ValueError("今天的最大晕质量必须大于最小晕质量。")
    if halo_track_point_count < 2:
        raise ValueError("halo_track_point_count 必须至少为 2。")

    model_name = concentration_model.lower()
    present_day_masses = np.geomspace(
        minimum_present_halo_mass_msun,
        maximum_present_halo_mass_msun,
        halo_track_point_count,
    )
    track_masses = np.asarray(
        mass_history_function(present_day_masses, z, model_name),
        dtype=float,
    )

    if model_name in {"prada12", "prada12_hmf_sigma"}:
        track_concentrations = np.asarray(
            hc.concentration_prada12_hmf_sigma(
                track_masses,
                z,
                cap_high_peak=True,
            ),
            dtype=float,
        )
    elif model_name == "ludlow16":
        track_concentrations = np.full(track_masses.shape, np.nan)
    else:
        raise ValueError(
            "concentration_model 必须是 'ludlow16'、'prada12' 或 "
            "'prada12_hmf_sigma'。"
        )

    rates_per_halo = np.asarray(
        [
            capture_rate_per_halo_a14_per_year(
                halo_mass_msun,
                z,
                mass_nodes_msun,
                mass_probability_weights,
                concentration_model=model_name,
                f_pbh=f_pbh,
                concentration_override=(
                    concentration
                    if model_name in {"prada12", "prada12_hmf_sigma"}
                    else None
                ),
            )
            for halo_mass_msun, concentration in zip(
                track_masses,
                track_concentrations,
            )
        ],
        dtype=float,
    )

    # 这是从原始 Fig. 8 反推出来的关键：HMF 网格均匀覆盖 M(z) 的端点，
    # 但它并不是逐点等于上面的非均匀 track_masses。
    paired_hmf_masses = np.geomspace(
        float(np.min(track_masses)),
        float(np.max(track_masses)),
        halo_track_point_count,
    )
    hmf_masses, dndlnm = press_schechter_halo_mass_function_hmfcalc(
        z,
        minimum_halo_mass_msun=float(np.min(track_masses)) * 0.8,
        maximum_halo_mass_msun=float(np.max(track_masses)) * 1.2,
        dlog10m=hmf_dlog10m,
    )
    finite_positive = (dndlnm > 0.0) & np.isfinite(dndlnm)
    paired_dndlnm = np.exp(
        np.interp(
            np.log(paired_hmf_masses),
            np.log(hmf_masses[finite_positive]),
            np.log(dndlnm[finite_positive]),
        )
    )

    rate_density_mpc3_per_year = np.trapezoid(
        rates_per_halo * paired_dndlnm,
        x=np.log(present_day_masses),
    )
    return float(rate_density_mpc3_per_year * MPC3_PER_GPC3)


def comoving_capture_rate_eq18_gpc3_per_year(
    z,
    mass_nodes_msun,
    mass_probability_weights,
    concentration_model="ludlow16",
    f_pbh=1.0,
    minimum_halo_mass_msun=1.0e3,
    maximum_halo_mass_msun=1.0e15,
    dlog10m=0.1,
):
    """计算论文 Eq. (18) 的总共动两体捕获率。

    论文写成

        R(z) = integral R_halo(M,z) * (dn/dM) dM.

    本函数用等价的对数质量形式计算：

        R(z) = integral R_halo(M,z) * (dn/dlnM) dlnM.

    ``R_halo`` 由 :func:`capture_rate_per_halo_a14_per_year` 计算，晕质量函数
    则由 :func:`press_schechter_halo_mass_function_hmfcalc` 调用 HMFcalc 的
    ``hmf`` 后端得到。最后把 ``Mpc^-3 yr^-1`` 乘以 ``10^9``，返回常用于
    引力波事件率的 ``Gpc^-3 yr^-1``。

    ``mass_nodes_msun`` 与 ``mass_probability_weights`` 的含义和 Appendix A14
    每晕率函数完全相同，所以单色质量函数传一个节点和权重 1，连续质量函数
    则传入数值积分节点和 ``psi(m)dm`` 权重。
    """
    halo_masses_msun, dndlnm_mpc3 = (
        press_schechter_halo_mass_function_hmfcalc(
            z,
            minimum_halo_mass_msun=minimum_halo_mass_msun,
            maximum_halo_mass_msun=maximum_halo_mass_msun,
            dlog10m=dlog10m,
        )
    )

    # 这里的每个 M 都表示红移 z 时实际存在的晕质量，直接代入 R_halo(M,z)。
    # 不再把它当作今天的 M0 后额外调用质量吸积史，否则会重复演化质量。
    rates_per_halo_per_year = np.asarray(
        [
            capture_rate_per_halo_a14_per_year(
                halo_mass_msun,
                z,
                mass_nodes_msun,
                mass_probability_weights,
                concentration_model=concentration_model,
                f_pbh=f_pbh,
            )
            for halo_mass_msun in halo_masses_msun
        ],
        dtype=float,
    )

    rate_density_mpc3_per_year = np.trapezoid(
        rates_per_halo_per_year * dndlnm_mpc3,
        x=np.log(halo_masses_msun),
    )
    return float(rate_density_mpc3_per_year * MPC3_PER_GPC3)
