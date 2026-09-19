"""Appendix D 中原初黑洞双星的初始轨道参数分布。

本模块实现目标论文 Appendix D 的 Eqs. (D1)--(D5)：

* 平均 PBH 间距 ``x_bar``；
* 重标度角动量的密度 ``P(j)``；
* 半长轴和角动量的联合密度 ``P(a, j)``；
* ``j = sqrt(1-e**2)`` 与偏心率 ``e`` 的相互转换；
* 在有限范围内采样相关的 ``(a, j, e)``。

公开长度统一使用 pc，质量统一使用 M_sun，密度统一使用 M_sun/pc^3。
论文只给出 ``z_eq ~= 3450``，没有在 Appendix D 单独列出计算 ``rho_eq``
所用的 ``Omega_m0`` 和 ``h``。因此本模块既允许调用者直接传入 ``rho_eq``，
也提供与项目 Ludlow16/Planck 参数一致的默认换算。

Fig. 19 的正文说逆变换采样使用 ``N_sample=10^4``，图注则写 ``N=10^5``。
本模块不替作者消除这个文字差异；样本数由调用者明确传入。

论文还没有说明“先采样 j”时 Eq. (D3) 中的 ``x/x_bar`` 取值，也没有说明
``0<j<1`` 与有限半长轴区间的归一化顺序。本模块采用数学上明确的做法：把
Eq. (D4) 作为联合密度，在指定矩形区域内整体截断，再保留 ``a`` 与 ``j`` 的
相关性进行采样。这一点在比较 Fig. 19 时必须与作者未公开的实现区别开来。
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid


# 与 haloconcentration.py 中 Ludlow16 的 Planck 宇宙学参数保持一致。
DEFAULT_OMEGA_M0 = 0.308
DEFAULT_H = 0.678
DEFAULT_Z_EQ = 3450.0
DEFAULT_SIGMA_EQ = 0.005
DEFAULT_ALPHA = 0.1

# G 的单位为 pc (km/s)^2 / M_sun；这样由 H0 和 G 算出的临界密度直接是
# M_sun/pc^3，不需要在函数中隐藏额外的 SI 单位转换。
GRAVITATIONAL_CONSTANT_PC_KM2_S2_MSUN = 4.30091727003628e-3
MEGAPARSEC_IN_PC = 1.0e6
AU_PER_PC = 206264.80624709636


def critical_density_z0_msun_pc3(h=DEFAULT_H):
    """返回今天的临界密度，单位为 ``M_sun/pc^3``。

    使用

    ``rho_crit,0 = 3 H0^2 / (8 pi G)``，

    其中 ``H0 = 100 h km s^-1 Mpc^-1``。
    """
    h = float(h)
    if h <= 0.0:
        raise ValueError("h 必须为正。")

    hubble_z0_km_s_pc = 100.0 * h / MEGAPARSEC_IN_PC
    return (
        3.0
        * hubble_z0_km_s_pc**2
        / (8.0 * np.pi * GRAVITATIONAL_CONSTANT_PC_KM2_S2_MSUN)
    )


def matter_density_at_equality_msun_pc3(
    z_eq=DEFAULT_Z_EQ,
    omega_m0=DEFAULT_OMEGA_M0,
    h=DEFAULT_H,
):
    """把今天的平均物质密度外推到物质-辐射相等时期。

    对非相对论物质使用 ``rho_m(z)=rho_m,0*(1+z)^3``。这里返回质量密度；
    Appendix D 在自然单位下把它写成平均能量密度 ``rho_eq``。
    """
    z_eq = float(z_eq)
    omega_m0 = float(omega_m0)
    if z_eq < 0.0:
        raise ValueError("z_eq 不能为负。")
    if omega_m0 <= 0.0:
        raise ValueError("omega_m0 必须为正。")

    rho_m0 = omega_m0 * critical_density_z0_msun_pc3(h)
    return rho_m0 * (1.0 + z_eq) ** 3


DEFAULT_RHO_EQ_MSUN_PC3 = matter_density_at_equality_msun_pc3()


def mean_pbh_separation_pc(
    pbh_mass_msun=30.0,
    f_pbh=1.0,
    rho_eq_msun_pc3=DEFAULT_RHO_EQ_MSUN_PC3,
):
    """计算 Appendix D Eq. (D1) 的平均 PBH 间距 ``x_bar``。

    ``x_bar = [3 m_PBH / (4 pi f_PBH rho_eq)]^(1/3)``。
    """
    pbh_mass_msun = np.asarray(pbh_mass_msun, dtype=float)
    f_pbh = float(f_pbh)
    rho_eq_msun_pc3 = float(rho_eq_msun_pc3)

    if np.any(pbh_mass_msun <= 0.0):
        raise ValueError("PBH 质量必须为正。")
    if f_pbh <= 0.0:
        raise ValueError("f_pbh 必须为正。")
    if rho_eq_msun_pc3 <= 0.0:
        raise ValueError("rho_eq_msun_pc3 必须为正。")

    return (
        3.0 * pbh_mass_msun / (4.0 * np.pi * f_pbh * rho_eq_msun_pc3)
    ) ** (1.0 / 3.0)


def separation_from_semimajor_axis_pc(
    semi_major_axis_pc,
    pbh_mass_msun=30.0,
    rho_eq_msun_pc3=DEFAULT_RHO_EQ_MSUN_PC3,
    alpha=DEFAULT_ALPHA,
):
    """计算 Appendix D Eq. (D5) 的初始 PBH 间距 ``x(a)``。

    ``x(a) = [3 a m_PBH / (4 pi alpha rho_eq)]^(1/4)``。
    """
    semi_major_axis_pc = np.asarray(semi_major_axis_pc, dtype=float)
    pbh_mass_msun = float(pbh_mass_msun)
    rho_eq_msun_pc3 = float(rho_eq_msun_pc3)
    alpha = float(alpha)

    if np.any(semi_major_axis_pc <= 0.0):
        raise ValueError("半长轴必须为正。")
    if pbh_mass_msun <= 0.0:
        raise ValueError("PBH 质量必须为正。")
    if rho_eq_msun_pc3 <= 0.0:
        raise ValueError("rho_eq_msun_pc3 必须为正。")
    if alpha <= 0.0:
        raise ValueError("alpha 必须为正。")

    return (
        3.0
        * semi_major_axis_pc
        * pbh_mass_msun
        / (4.0 * np.pi * alpha * rho_eq_msun_pc3)
    ) ** 0.25


def angular_momentum_scale(
    x_over_xbar,
    f_pbh=1.0,
    sigma_eq=DEFAULT_SIGMA_EQ,
):
    """返回 Eq. (D3) 分母中的角动量尺度 ``j_0``。

    写成 ``y(j)=j/j_0`` 后，

    ``j_0 = 0.5*sqrt(1 + sigma_eq^2/f_pbh^2)*(x/x_bar)^3``。
    """
    x_over_xbar = np.asarray(x_over_xbar, dtype=float)
    f_pbh = float(f_pbh)
    sigma_eq = float(sigma_eq)

    if np.any(x_over_xbar <= 0.0):
        raise ValueError("x_over_xbar 必须为正。")
    if f_pbh <= 0.0:
        raise ValueError("f_pbh 必须为正。")
    if sigma_eq < 0.0:
        raise ValueError("sigma_eq 不能为负。")

    return (
        0.5
        * np.sqrt(1.0 + sigma_eq**2 / f_pbh**2)
        * x_over_xbar**3
    )


def angular_momentum_cdf(j, j_scale):
    """返回 Eq. (D2) 从 0 到 ``j`` 的解析积分。

    若 ``y=j/j_scale``，则

    ``integral_0^j P(j') dj' = 1 - 1/sqrt(1+y^2)``。

    这个结果在 ``j -> infinity`` 时趋于 1。论文采样限制在 ``0<j<1``，
    因此有限区间采样时还要除以 ``CDF(1)``。
    """
    j = np.asarray(j, dtype=float)
    j_scale = np.asarray(j_scale, dtype=float)
    if np.any(j < 0.0):
        raise ValueError("j 不能为负。")
    if np.any(j_scale <= 0.0):
        raise ValueError("j_scale 必须为正。")

    y = j / j_scale
    return 1.0 - 1.0 / np.sqrt(1.0 + y**2)


def angular_momentum_pdf(
    j,
    x_over_xbar,
    f_pbh=1.0,
    sigma_eq=DEFAULT_SIGMA_EQ,
    normalize_on_unit_interval=False,
):
    """计算 Appendix D Eqs. (D2)--(D3) 的 ``P(j)``。

    默认返回论文印刷的原始密度。设置
    ``normalize_on_unit_interval=True`` 时，返回在 ``0<j<1`` 上重新归一化的
    条件密度，适合直接与有限区间直方图比较。
    """
    j = np.asarray(j, dtype=float)
    j_scale = angular_momentum_scale(x_over_xbar, f_pbh, sigma_eq)
    j, j_scale = np.broadcast_arrays(j, j_scale)

    density = np.zeros(j.shape, dtype=float)
    inside = (j > 0.0) & (j < 1.0)
    if np.any(inside):
        y = j[inside] / j_scale[inside]
        density[inside] = (
            y**2 / (j[inside] * (1.0 + y**2) ** 1.5)
        )

    if normalize_on_unit_interval:
        unit_interval_probability = angular_momentum_cdf(1.0, j_scale)
        density = density / unit_interval_probability

    return density


def joint_orbital_pdf(
    semi_major_axis_pc,
    j,
    pbh_mass_msun=30.0,
    f_pbh=1.0,
    rho_eq_msun_pc3=DEFAULT_RHO_EQ_MSUN_PC3,
    sigma_eq=DEFAULT_SIGMA_EQ,
    alpha=DEFAULT_ALPHA,
):
    """计算 Appendix D Eq. (D4) 的原始联合密度 ``P(a,j)``。

    返回值对 ``a`` 的单位是 ``pc^-1``。这个函数忠实返回印刷公式，没有在
    Fig. 19 的有限 ``a``、``j`` 范围上另做整体归一化；采样器会单独处理截断。
    """
    semi_major_axis_pc = np.asarray(semi_major_axis_pc, dtype=float)
    j = np.asarray(j, dtype=float)
    semi_major_axis_pc, j = np.broadcast_arrays(semi_major_axis_pc, j)

    if np.any(semi_major_axis_pc <= 0.0):
        raise ValueError("半长轴必须为正。")

    x_pc = separation_from_semimajor_axis_pc(
        semi_major_axis_pc,
        pbh_mass_msun,
        rho_eq_msun_pc3,
        alpha,
    )
    x_bar_pc = mean_pbh_separation_pc(
        pbh_mass_msun,
        f_pbh,
        rho_eq_msun_pc3,
    )
    x_over_xbar = x_pc / x_bar_pc

    semi_major_axis_factor = (
        3.0
        * semi_major_axis_pc ** (-0.25)
        / 4.0
        * (f_pbh / (alpha * x_bar_pc)) ** 0.75
        * np.exp(-(x_over_xbar**3))
    )
    return semi_major_axis_factor * angular_momentum_pdf(
        j,
        x_over_xbar,
        f_pbh=f_pbh,
        sigma_eq=sigma_eq,
        normalize_on_unit_interval=False,
    )


def angular_momentum_from_eccentricity(eccentricity):
    """由 ``e`` 计算 ``j=sqrt(1-e^2)``。"""
    eccentricity = np.asarray(eccentricity, dtype=float)
    if np.any((eccentricity < 0.0) | (eccentricity >= 1.0)):
        raise ValueError("偏心率必须满足 0 <= e < 1。")
    return np.sqrt(1.0 - eccentricity**2)


def eccentricity_from_angular_momentum(j):
    """由 ``j`` 计算 ``e=sqrt(1-j^2)``。"""
    j = np.asarray(j, dtype=float)
    if np.any((j <= 0.0) | (j >= 1.0)):
        raise ValueError("角动量必须满足 0 < j < 1。")
    return np.sqrt(1.0 - j**2)


def semi_major_axis_marginal_pdf(
    semi_major_axis_pc,
    pbh_mass_msun=30.0,
    f_pbh=1.0,
    rho_eq_msun_pc3=DEFAULT_RHO_EQ_MSUN_PC3,
    sigma_eq=DEFAULT_SIGMA_EQ,
    alpha=DEFAULT_ALPHA,
):
    """把 Eq. (D4) 在 ``0<j<1`` 上积分，得到 ``a`` 的未归一化边缘密度。"""
    semi_major_axis_pc = np.asarray(semi_major_axis_pc, dtype=float)
    if np.any(semi_major_axis_pc <= 0.0):
        raise ValueError("半长轴必须为正。")

    x_pc = separation_from_semimajor_axis_pc(
        semi_major_axis_pc,
        pbh_mass_msun,
        rho_eq_msun_pc3,
        alpha,
    )
    x_bar_pc = mean_pbh_separation_pc(
        pbh_mass_msun,
        f_pbh,
        rho_eq_msun_pc3,
    )
    x_over_xbar = x_pc / x_bar_pc
    j_scale = angular_momentum_scale(x_over_xbar, f_pbh, sigma_eq)

    semi_major_axis_factor = (
        3.0
        * semi_major_axis_pc ** (-0.25)
        / 4.0
        * (f_pbh / (alpha * x_bar_pc)) ** 0.75
        * np.exp(-(x_over_xbar**3))
    )
    return semi_major_axis_factor * angular_momentum_cdf(1.0, j_scale)


@dataclass(frozen=True)
class OrbitalSamples:
    """一次采样得到的三个相互对应的轨道参数数组。"""

    semi_major_axis_pc: np.ndarray
    angular_momentum: np.ndarray
    eccentricity: np.ndarray


def sample_initial_orbital_parameters(
    sample_count,
    rng,
    pbh_mass_msun=30.0,
    f_pbh=1.0,
    rho_eq_msun_pc3=DEFAULT_RHO_EQ_MSUN_PC3,
    sigma_eq=DEFAULT_SIGMA_EQ,
    alpha=DEFAULT_ALPHA,
    minimum_semi_major_axis_pc=1.0e-6,
    maximum_semi_major_axis_pc=1.0,
    semi_major_axis_grid_point_count=16384,
):
    """从截断后的 Appendix D 联合分布采样 ``(a,j,e)``。

    算法分两步：

    1. 在对数 ``a`` 网格上计算 ``P(a,j)`` 对 ``0<j<1`` 的解析积分，随后对
       ``a`` 的数值 CDF 做逆变换采样；
    2. 对每个已经抽到的 ``a``，使用 Eq. (D2) 的解析 CDF 在 ``0<j<1`` 上
       采样条件角动量，再计算 ``e=sqrt(1-j^2)``。

    这样保留了 ``a`` 与 ``j`` 的相关性。``rng`` 必须由调用者显式提供，例如
    ``np.random.default_rng(12345)``；本函数不会修改 NumPy 的全局随机状态。
    """
    sample_count = int(sample_count)
    semi_major_axis_grid_point_count = int(semi_major_axis_grid_point_count)
    minimum_semi_major_axis_pc = float(minimum_semi_major_axis_pc)
    maximum_semi_major_axis_pc = float(maximum_semi_major_axis_pc)

    if sample_count <= 0:
        raise ValueError("sample_count 必须为正整数。")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng 必须是 numpy.random.Generator。")
    if minimum_semi_major_axis_pc <= 0.0:
        raise ValueError("最小半长轴必须为正。")
    if maximum_semi_major_axis_pc <= minimum_semi_major_axis_pc:
        raise ValueError("最大半长轴必须大于最小半长轴。")
    if semi_major_axis_grid_point_count < 2:
        raise ValueError("半长轴 CDF 网格至少需要两个点。")

    semi_major_axis_grid_pc = np.geomspace(
        minimum_semi_major_axis_pc,
        maximum_semi_major_axis_pc,
        semi_major_axis_grid_point_count,
    )
    marginal_density = semi_major_axis_marginal_pdf(
        semi_major_axis_grid_pc,
        pbh_mass_msun,
        f_pbh,
        rho_eq_msun_pc3,
        sigma_eq,
        alpha,
    )
    semi_major_axis_cdf = cumulative_trapezoid(
        marginal_density,
        x=semi_major_axis_grid_pc,
        initial=0.0,
    )
    semi_major_axis_cdf /= semi_major_axis_cdf[-1]

    semi_major_axis_pc = np.interp(
        rng.random(sample_count),
        semi_major_axis_cdf,
        semi_major_axis_grid_pc,
    )

    x_pc = separation_from_semimajor_axis_pc(
        semi_major_axis_pc,
        pbh_mass_msun,
        rho_eq_msun_pc3,
        alpha,
    )
    x_bar_pc = mean_pbh_separation_pc(
        pbh_mass_msun,
        f_pbh,
        rho_eq_msun_pc3,
    )
    j_scale = angular_momentum_scale(x_pc / x_bar_pc, f_pbh, sigma_eq)

    # F(1) 是未截断 D2 分布落在 0<j<1 内的概率。先将 [0,1) 上的均匀
    # 随机数缩放到 [0,F(1))，再使用解析逆 CDF。
    probability_below_one = angular_momentum_cdf(1.0, j_scale)
    target_cdf = rng.random(sample_count) * probability_below_one
    inverse_cdf_denominator = 1.0 - target_cdf
    y = np.sqrt(inverse_cdf_denominator ** (-2.0) - 1.0)
    angular_momentum = np.minimum(j_scale * y, np.nextafter(1.0, 0.0))
    eccentricity = eccentricity_from_angular_momentum(angular_momentum)

    return OrbitalSamples(
        semi_major_axis_pc=semi_major_axis_pc,
        angular_momentum=angular_momentum,
        eccentricity=eccentricity,
    )


@dataclass(frozen=True)
class AppendixDOrbitalDistribution:
    """保存 Appendix D 参数，并提供计划约定的 ``sample(n, rng)`` 接口。"""

    pbh_mass_msun: float = 30.0
    f_pbh: float = 1.0
    rho_eq_msun_pc3: float = DEFAULT_RHO_EQ_MSUN_PC3
    sigma_eq: float = DEFAULT_SIGMA_EQ
    alpha: float = DEFAULT_ALPHA
    minimum_semi_major_axis_pc: float = 1.0e-6
    maximum_semi_major_axis_pc: float = 1.0
    semi_major_axis_grid_point_count: int = 16384

    def sample(self, sample_count, rng):
        """调用有限区间联合采样器，并返回 ``OrbitalSamples``。"""
        return sample_initial_orbital_parameters(
            sample_count=sample_count,
            rng=rng,
            pbh_mass_msun=self.pbh_mass_msun,
            f_pbh=self.f_pbh,
            rho_eq_msun_pc3=self.rho_eq_msun_pc3,
            sigma_eq=self.sigma_eq,
            alpha=self.alpha,
            minimum_semi_major_axis_pc=self.minimum_semi_major_axis_pc,
            maximum_semi_major_axis_pc=self.maximum_semi_major_axis_pc,
            semi_major_axis_grid_point_count=(
                self.semi_major_axis_grid_point_count
            ),
        )
