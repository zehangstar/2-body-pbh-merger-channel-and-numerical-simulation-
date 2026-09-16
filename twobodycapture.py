"""暗物质晕中的 PBH 两体引力波捕获。

本模块实现目标论文 Eq. (11) 的捕获截面、速度平均 ``<Sigma v>``，以及
Appendix A14 对任意 PBH 质量分布的每晕捕获率。单色和连续质量谱使用同一个
每晕率函数，区别只在传入的质量积分节点和权重。

速度分布沿用 ``halostructure.py`` 中归一化的三维各向同性密度 ``f_3(v)``：

    4*pi*integral v^2*f_3(v) dv = 1.

因此速度平均显式使用球壳测度 ``4*pi*v^2 dv``，不会重复或遗漏 Jacobian。
这里没有自动自检、断言或独立测试模块。
"""

import numpy as np
from scipy.integrate import quad

import haloconcentration as hc
import halostructure as hs


SPEED_OF_LIGHT_KM_S = 299792.458
SECONDS_PER_YEAR = 365.25 * 24.0 * 3600.0
KILOMETERS_PER_KPC = 3.0856775814913673e16
KM_S_TO_KPC_PER_YEAR = SECONDS_PER_YEAR / KILOMETERS_PER_KPC


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
):
    """用 Appendix A14 计算任意 PBH 质量分布的每晕捕获率，单位 yr^-1。

    ``mass_nodes_msun`` 与 ``mass_probability_weights`` 共同表示质量测度：

        integral psi(m) F(m) dm ~= sum_i weight_i F(m_i).

    单色谱传入一个节点和权重 1；连续谱传入数值积分节点与 ``psi(m)dm``
    权重。两者调用完全相同的本函数，不存在单色专用每晕率公式。

    ``f_pbh`` 是参与直接捕获的有效 PBH 质量占暗物质的比例。论文基准在
    Eq. (12) 附近取总 PBH 丰度为 1；若只允许部分 PBH 参与，应显式传入相应值，
    结果按 ``f_pbh**2`` 缩放。
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
    elif model_name == "prada12":
        cosmology = hc.PRADA_COSMOLOGY
    else:
        raise ValueError("concentration_model 必须是 'ludlow16' 或 'prada12'。")

    concentration = float(
        hc.concentration_model(
            halo_mass_msun,
            z,
            model_name,
            cap_prada=(model_name == "prada12"),
        )
    )
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
