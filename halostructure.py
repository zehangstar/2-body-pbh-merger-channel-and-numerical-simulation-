"""NFW 暗物质晕结构：9 月 10 日阶段的逐步实现。

本模块接收前两阶段已经得到的晕质量、浓度和红移，逐步构造：

    (M, C, z) -> (R_vir, R_s, rho_s) -> rho_NFW(r) -> M(<r)

当前模块实现 NFW 特征密度、密度剖面、内部质量、特征速度和论文采用的
截断 Maxwell 速率分布；两体捕获率将在后续阶段加入。

这里没有自动自检、断言或独立测试模块。
"""

import numpy as np
from scipy.special import erf

import haloconcentration as hc


GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN = 4.30091e-6


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


def truncated_maxwell_speed_pdf(
    speed_km_s, velocity_dispersion_km_s, cutoff_velocity_km_s
):
    """返回归一化的 lowered-Maxwell 速率 PDF，单位为 (km/s)^-1。

    分布在 ``0 <= speed <= cutoff_velocity`` 上正比于
    ``speed**2 * (exp(-speed**2 / dispersion**2)
    - exp(-cutoff**2 / dispersion**2))``，在区间外为零。
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
    normalization = velocity_dispersion_km_s**3 * (
        0.25 * np.sqrt(np.pi) * erf(q)
        - np.exp(-(q**2)) * (0.5 * q + q**3 / 3.0)
    )
    profile = speed_km_s**2 * (
        np.exp(-(speed_km_s / velocity_dispersion_km_s) ** 2)
        - np.exp(-(q**2))
    )
    return np.where(
        (speed_km_s >= 0.0) & (speed_km_s <= cutoff_velocity_km_s),
        profile / normalization,
        0.0,
    )
