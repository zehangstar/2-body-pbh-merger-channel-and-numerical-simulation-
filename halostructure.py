"""NFW 暗物质晕结构：9 月 10 日阶段的逐步实现。

本模块接收前两阶段已经得到的晕质量、浓度和红移，逐步构造：

    (M, C, z) -> (R_vir, R_s, rho_s) -> rho_NFW(r) -> M(<r)

当前模块实现 NFW 特征密度、密度剖面、内部质量、特征速度和 Bird et al.
采用的截断 Maxwell 三维速度密度；并按照论文 Eqs. (23)、(24)、(27)、(28)
构造 binary-single 演化所需的动态晕壳层环境。

这里没有自动自检、断言或独立测试模块。
"""

from dataclasses import dataclass

import numpy as np
from scipy.special import erf

import haloconcentration as hc


GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN = 4.30091e-6
KILOPARSEC_IN_PC = 1000.0
AU_PER_PC = 206264.80624709636
VIRIAL_OVERDENSITY = 200.0


@dataclass(frozen=True)
class HaloShellState:
    """一个现今质量标签 ``M0`` 在红移 ``z`` 的完整壳层环境。

    数组下标始终代表同一个球壳。该对象只提供 Eqs. (19)、(25) 后续会消费的
    环境量，不在这里推进双星半长轴、偏心率或记录并合事件。
    """

    present_day_halo_mass_msun: float
    redshift: float
    halo_mass_msun: float
    concentration: float
    virial_radius_pc: float
    scale_radius_pc: float
    scale_density_msun_pc3: float
    shell_boundaries_pc: np.ndarray
    shell_midpoints_pc: np.ndarray
    velocity_evaluation_radii_pc: np.ndarray
    shell_density_msun_pc3: np.ndarray
    enclosed_mass_at_velocity_radius_msun: np.ndarray
    shell_mass_msun: np.ndarray
    velocity_dispersion_km_s: np.ndarray
    hard_semimajor_axis_pc: np.ndarray
    hard_semimajor_axis_au: np.ndarray


def nfw_characteristic_overdensity(concentration, virial_overdensity=200.0):
    """返回 NFW 尺度密度相对于临界密度的无量纲比值。

    论文 Eqs. (2)--(3) 为

        delta_NFW(C) = (Delta / 3) * C**3 / g(C),
        g(C) = ln(1 + C) - C/(1 + C),

    其中论文采用 ``Delta=200``。因此尺度密度最终可以写成

        rho_s(z) = rho_crit(z) * delta_NFW(C).

    Parameters
    ----------
    concentration : float 或 array-like
        NFW 浓度 C=R_vir/R_s。它是无量纲量。
    virial_overdensity : float, optional
        virial 半径定义中的平均密度倍数 Delta，默认使用论文的 200。

    Returns
    -------
    delta_nfw : numpy.ndarray
        NFW 特征密度对比，是无量纲量。输入可以是单个浓度或浓度数组。

    Notes
    -----
    这里的 ``delta_nfw`` 不是峰高公式中的球形坍缩阈值
    ``delta_sc = 1.686``。两者虽然都常被文献写成 delta_c，却描述完全不同的量。
    """
    concentration = np.asarray(concentration, dtype=float)

    # g(C) 已在 haloconcentration.py 中实现。直接调用它可以确保浓度模块、
    # 质量吸积模块与这里使用完全相同的 NFW 质量积分函数，而不是复制一份公式。
    mass_integral_factor = hc.nfw_g(concentration)

    # 从 M_vir=4*pi*rho_s*R_s^3*g(C) 与
    # M_vir=(4*pi/3)*Delta*rho_crit*(C*R_s)^3 对比即可得到本式。
    return (
        virial_overdensity
        * concentration**3
        / (3.0 * mass_integral_factor)
    )


def rho_crit(z, cosmology=hc.LUDLOW_COSMOLOGY):
    """返回红移 z 处的宇宙临界密度，单位为 M_sun/kpc^3。"""
    z = np.asarray(z, dtype=float)
    hubble0_km_s_kpc = 0.1 * cosmology["h"]
    expansion_rate_squared = (
        cosmology["omega_m0"] * (1.0 + z) ** 3
        + cosmology["omega_lambda0"]
    )
    return (
        3.0
        * hubble0_km_s_kpc**2
        * expansion_rate_squared
        / (8.0 * np.pi * GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN)
    )


def rho_s(
    concentration,
    z,
    virial_overdensity=200.0,
    cosmology=hc.LUDLOW_COSMOLOGY,
):
    """返回 NFW 尺度密度 rho_s(C,z)，单位为 M_sun/kpc^3。"""
    return rho_crit(z, cosmology) * nfw_characteristic_overdensity(
        concentration, virial_overdensity
    )


def rho_nfw(radius_kpc, scale_radius_kpc, scale_density_msun_kpc3):
    """返回半径 ``radius_kpc`` 处的 NFW 密度，单位为 M_sun/kpc^3。"""
    radius_kpc = np.asarray(radius_kpc, dtype=float)
    scale_radius_kpc = np.asarray(scale_radius_kpc, dtype=float)
    scale_density_msun_kpc3 = np.asarray(
        scale_density_msun_kpc3, dtype=float
    )
    if np.any(radius_kpc < 0.0) or np.any(scale_radius_kpc <= 0.0):
        raise ValueError("radius 必须非负，scale radius 必须为正。")

    x = radius_kpc / scale_radius_kpc
    with np.errstate(divide="ignore"):
        return scale_density_msun_kpc3 / (x * (1.0 + x) ** 2)


def nfw_enclosed_mass(
    radius_kpc, scale_radius_kpc, scale_density_msun_kpc3
):
    """返回 NFW 晕在 ``radius_kpc`` 内的质量，单位为 M_sun。"""
    radius_kpc = np.asarray(radius_kpc, dtype=float)
    scale_radius_kpc = np.asarray(scale_radius_kpc, dtype=float)
    scale_density_msun_kpc3 = np.asarray(
        scale_density_msun_kpc3, dtype=float
    )
    if np.any(radius_kpc < 0.0) or np.any(scale_radius_kpc <= 0.0):
        raise ValueError("radius 必须非负，scale radius 必须为正。")

    x = radius_kpc / scale_radius_kpc
    return (
        4.0
        * np.pi
        * scale_density_msun_kpc3
        * scale_radius_kpc**3
        * hc.nfw_g(x)
    )


def nfw_characteristic_velocity(
    scale_radius_kpc,
    scale_density_msun_kpc3,
    x_max=2.1626,
):
    """返回 NFW 最大圆周速度 ``V_max``，单位为 km/s。"""
    scale_radius_kpc = np.asarray(scale_radius_kpc, dtype=float)
    scale_density_msun_kpc3 = np.asarray(
        scale_density_msun_kpc3, dtype=float
    )
    if np.any(scale_radius_kpc <= 0.0) or x_max <= 0.0:
        raise ValueError("scale radius 和 x_max 必须为正。")

    return np.sqrt(
        4.0
        * np.pi
        * GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN
        * scale_density_msun_kpc3
        * scale_radius_kpc**2
        * hc.nfw_g(x_max)
        / x_max
    )


def virial_cutoff_velocity(halo_mass_msun, virial_radius_kpc):
    """返回论文所用截止速度 sqrt(2 G M_vir / R_vir)，单位为 km/s。"""
    halo_mass_msun = np.asarray(halo_mass_msun, dtype=float)
    virial_radius_kpc = np.asarray(virial_radius_kpc, dtype=float)
    if np.any(halo_mass_msun < 0.0) or np.any(virial_radius_kpc <= 0.0):
        raise ValueError("halo mass 必须非负，virial radius 必须为正。")

    return np.sqrt(
        2.0
        * GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN
        * halo_mass_msun
        / virial_radius_kpc
    )


def paper_shell_count(present_day_halo_mass_msun):
    """返回论文 Section IV.B 对现今晕质量采用的球壳数。

    论文按现今质量的数量级设置 1、2、3、5、10 个壳层。壳层数沿同一条
    ``M0`` 吸积轨迹保持不变；不能随 ``M(z)`` 跨过阈值而重新分壳。
    """
    present_day_halo_mass_msun = float(present_day_halo_mass_msun)
    if present_day_halo_mass_msun <= 0.0:
        raise ValueError("现今晕质量必须为正。")
    if present_day_halo_mass_msun < 1.0e4:
        return 1
    if present_day_halo_mass_msun < 1.0e5:
        return 2
    if present_day_halo_mass_msun < 1.0e6:
        return 3
    if present_day_halo_mass_msun < 1.0e7:
        return 5
    return 10


def logarithmic_shell_boundaries_pc(virial_radius_pc, shell_count):
    """按论文 Eq. (27) 返回从 0 到 ``R_vir`` 的壳边界，单位 pc。"""
    virial_radius_pc = float(virial_radius_pc)
    shell_count = int(shell_count)
    if virial_radius_pc <= 0.0:
        raise ValueError("virial radius 必须为正。")
    if shell_count <= 0:
        raise ValueError("shell_count 必须为正整数。")

    boundary_indices = np.arange(shell_count + 1, dtype=float)
    return np.expm1(
        boundary_indices
        / shell_count
        * np.log1p(virial_radius_pc)
    )


def binary_single_velocity_dispersion_km_s(
    radius_kpc,
    enclosed_mass_msun,
):
    """计算论文 Eq. (24) 的局部环境速度弥散 ``sqrt(2GM(<r)/r)``。

    原文 Table III 的 ``a_h`` 数值更接近把这里的速度再除以 2 的约定；
    本函数按印刷 Eq. (24) 实现，不用表格数值反向改写公式。
    """
    radius_kpc = np.asarray(radius_kpc, dtype=float)
    enclosed_mass_msun = np.asarray(enclosed_mass_msun, dtype=float)
    if np.any(radius_kpc <= 0.0):
        raise ValueError("计算局部速度弥散的半径必须为正。")
    if np.any(enclosed_mass_msun < 0.0):
        raise ValueError("包围质量不能为负。")
    return np.sqrt(
        2.0
        * GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN
        * enclosed_mass_msun
        / radius_kpc
    )


def hard_binary_semimajor_axis_pc(
    primary_pbh_mass_msun,
    velocity_dispersion_km_s,
):
    """计算论文 Eq. (23) 的硬双星临界半长轴，返回 pc。"""
    primary_pbh_mass_msun = float(primary_pbh_mass_msun)
    velocity_dispersion_km_s = np.asarray(
        velocity_dispersion_km_s,
        dtype=float,
    )
    if primary_pbh_mass_msun <= 0.0:
        raise ValueError("primary PBH mass 必须为正。")
    if np.any(velocity_dispersion_km_s <= 0.0):
        raise ValueError("局部速度弥散必须为正。")
    return (
        GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN
        * primary_pbh_mass_msun
        / (4.0 * velocity_dispersion_km_s**2)
        * KILOPARSEC_IN_PC
    )


def build_halo_shell_state(
    present_day_halo_mass_msun,
    z,
    primary_pbh_mass_msun=30.0,
    model=hc.DEFAULT_HALO_MODEL,
    shell_count=None,
):
    """构造一个 ``(M0,z)`` 对应的 binary-single 晕壳层环境。

    默认计算链为

    ``M0 -> Prada12-HMF M(z) -> C[M(z),z] -> NFW -> shell state``。

    ``shell_count`` 默认按论文的现今质量分级选择，也可为诊断显式覆盖。
    密度按 Eq. (28) 在壳中点取值。Eq. (24) 按原文 Table III 的位置口径：
    单球模型取 ``R_vir/2``，多壳模型在各壳外边界 ``R_i`` 求
    ``M(<R_i)`` 和 ``v_disp``。两个评价半径都显式保存在返回对象中，后续
    Eqs. (19)、(25) 不需要猜测数组含义。
    """
    present_day_halo_mass_msun = float(present_day_halo_mass_msun)
    z = float(z)
    if present_day_halo_mass_msun <= 0.0:
        raise ValueError("现今晕质量必须为正。")
    if z < 0.0:
        raise ValueError("红移 z 不能为负。")

    if shell_count is None:
        shell_count = paper_shell_count(present_day_halo_mass_msun)
    else:
        shell_count = int(shell_count)
        if shell_count <= 0:
            raise ValueError("shell_count 必须为正整数。")

    cosmology = hc.cosmology_for_model(model)
    halo_mass_msun = float(
        hc.mass_accretion_history(
            present_day_halo_mass_msun,
            z,
            model=model,
        )
    )
    concentration = float(
        hc.concentration_model(
            halo_mass_msun,
            z,
            model=model,
        )
    )
    critical_density_msun_kpc3 = float(rho_crit(z, cosmology))
    virial_radius_kpc = (
        3.0
        * halo_mass_msun
        / (
            4.0
            * np.pi
            * VIRIAL_OVERDENSITY
            * critical_density_msun_kpc3
        )
    ) ** (1.0 / 3.0)
    scale_radius_kpc = virial_radius_kpc / concentration
    scale_density_msun_kpc3 = float(
        rho_s(concentration, z, cosmology=cosmology)
    )

    virial_radius_pc = virial_radius_kpc * KILOPARSEC_IN_PC
    boundaries_pc = logarithmic_shell_boundaries_pc(
        virial_radius_pc,
        shell_count,
    )
    midpoints_pc = 0.5 * (boundaries_pc[:-1] + boundaries_pc[1:])
    boundaries_kpc = boundaries_pc / KILOPARSEC_IN_PC
    midpoints_kpc = midpoints_pc / KILOPARSEC_IN_PC
    velocity_radii_pc = (
        midpoints_pc
        if shell_count == 1
        else boundaries_pc[1:]
    )
    velocity_radii_kpc = velocity_radii_pc / KILOPARSEC_IN_PC

    shell_density_msun_kpc3 = rho_nfw(
        midpoints_kpc,
        scale_radius_kpc,
        scale_density_msun_kpc3,
    )
    enclosed_at_velocity_radius_msun = nfw_enclosed_mass(
        velocity_radii_kpc,
        scale_radius_kpc,
        scale_density_msun_kpc3,
    )
    enclosed_at_boundaries_msun = nfw_enclosed_mass(
        boundaries_kpc,
        scale_radius_kpc,
        scale_density_msun_kpc3,
    )
    shell_mass_msun = np.diff(enclosed_at_boundaries_msun)
    velocity_dispersion_km_s = binary_single_velocity_dispersion_km_s(
        velocity_radii_kpc,
        enclosed_at_velocity_radius_msun,
    )
    hard_semimajor_axis_pc = hard_binary_semimajor_axis_pc(
        primary_pbh_mass_msun,
        velocity_dispersion_km_s,
    )

    return HaloShellState(
        present_day_halo_mass_msun=present_day_halo_mass_msun,
        redshift=z,
        halo_mass_msun=halo_mass_msun,
        concentration=concentration,
        virial_radius_pc=virial_radius_pc,
        scale_radius_pc=scale_radius_kpc * KILOPARSEC_IN_PC,
        scale_density_msun_pc3=(
            scale_density_msun_kpc3 / KILOPARSEC_IN_PC**3
        ),
        shell_boundaries_pc=boundaries_pc,
        shell_midpoints_pc=midpoints_pc,
        velocity_evaluation_radii_pc=velocity_radii_pc,
        shell_density_msun_pc3=(
            shell_density_msun_kpc3 / KILOPARSEC_IN_PC**3
        ),
        enclosed_mass_at_velocity_radius_msun=(
            enclosed_at_velocity_radius_msun
        ),
        shell_mass_msun=shell_mass_msun,
        velocity_dispersion_km_s=velocity_dispersion_km_s,
        hard_semimajor_axis_pc=hard_semimajor_axis_pc,
        hard_semimajor_axis_au=hard_semimajor_axis_pc * AU_PER_PC,
    )


def build_halo_shell_history(
    present_day_halo_mass_msun,
    redshifts,
    primary_pbh_mass_msun=30.0,
    model=hc.DEFAULT_HALO_MODEL,
    shell_count=None,
):
    """返回同一 ``M0`` 轨迹在多个红移上的 ``HaloShellState`` 元组。"""
    redshifts = np.atleast_1d(
        np.asarray(redshifts, dtype=float)
    ).ravel()
    if np.any(redshifts < 0.0):
        raise ValueError("红移 z 不能为负。")

    fixed_shell_count = (
        paper_shell_count(present_day_halo_mass_msun)
        if shell_count is None
        else int(shell_count)
    )
    return tuple(
        build_halo_shell_state(
            present_day_halo_mass_msun=present_day_halo_mass_msun,
            z=redshift,
            primary_pbh_mass_msun=primary_pbh_mass_msun,
            model=model,
            shell_count=fixed_shell_count,
        )
        for redshift in redshifts
    )


def truncated_maxwell_velocity_density(
    speed_km_s, velocity_dispersion_km_s, cutoff_velocity_km_s
):
    """返回 Bird et al. Eq. (7) 的三维各向同性速度密度。

    返回值的单位为 ``(km/s)^-3``，并满足
    ``4*pi*integral(P(v)*v**2*dv) = 1``。因此这里的 ``P(v)`` 本身
    不含速度空间球壳因子 ``v**2``；该因子应在后续速度积分的测度中加入。
    """
    speed_km_s = np.asarray(speed_km_s, dtype=float)
    velocity_dispersion_km_s = np.asarray(
        velocity_dispersion_km_s, dtype=float
    )
    cutoff_velocity_km_s = np.asarray(cutoff_velocity_km_s, dtype=float)
    if np.any(velocity_dispersion_km_s <= 0.0) or np.any(
        cutoff_velocity_km_s <= 0.0
    ):
        raise ValueError("velocity dispersion 和 cutoff velocity 必须为正。")

    q = cutoff_velocity_km_s / velocity_dispersion_km_s
    radial_integral = velocity_dispersion_km_s**3 * (
        0.25 * np.sqrt(np.pi) * erf(q)
        - np.exp(-(q**2)) * (0.5 * q + q**3 / 3.0)
    )
    normalization_coefficient = 1.0 / (4.0 * np.pi * radial_integral)
    profile = (
        np.exp(-(speed_km_s / velocity_dispersion_km_s) ** 2)
        - np.exp(-(q**2))
    )
    return np.where(
        (speed_km_s >= 0.0) & (speed_km_s <= cutoff_velocity_km_s),
        normalization_coefficient * profile,
        0.0,
    )
