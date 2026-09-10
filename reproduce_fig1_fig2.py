"""复现 Aljaf & Cholis (2025) 的 Fig. 1 和 Fig. 2。

学习约定
--------
这个文件最初采用逐函数填空的学习方式。为了在 2026-09-20 前完成论文的
全流程精读与分层复现，从本版本开始改为“模块交付 + 逐段精讲”：函数保持
短小、注释解释公式和单位，同时每个阶段直接形成能够生成论文图像的程序。

当前版本包括：
1. Ludlow16 Appendix C 的浓度拟合；
2. Prada12 Eqs. (12)--(23) 的浓度拟合；
3. PBH 论文 Appendix C 的平均质量吸积史；
4. Fig. 1 和 Fig. 2 的绘制及曲线数据导出。

公开接口中的 halo mass 一律使用物理太阳质量 M_sun。原始拟合公式若用
h^-1 M_sun，会在函数内部显式转换，避免单位被悄悄混用。

给 Python 初学者的阅读提示
-------------------------
1. ``def function_name(...):`` 表示定义函数，括号内是函数的输入。
2. ``return result`` 表示把计算结果交还给调用函数的位置。没有 return 的函数
   默认返回 None，后面的物理计算便无法继续。
3. Python 用 ``**`` 表示乘方。例如 ``x**3`` 是 x 的三次方，不要写 ``x^3``；
   在 Python 中 ``^`` 是按位异或，不是数学乘方。
4. NumPy 数组可以一次计算许多红移。例如 ``1.0 + z`` 在 z 是数组时会对每个
   元素分别加 1，因此这里不需要写 for 循环。
5. 所有公开质量输入均为物理质量 M_sun；变量名带 ``_msun`` 是单位提醒。
6. 本文件不设置独立的自检函数或 assert；数值结果通过输出表和论文图像人工核对。
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib

# 当前环境没有可用的 Tk 图形界面。Agg 后端直接把图写入 PNG，不弹出窗口，
# 也不会改变任何物理计算。
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import quad


# Ludlow16 Appendix C 所用 Planck cosmology。
LUDLOW_COSMOLOGY = {
    "name": "Planck (Ludlow16)",
    "omega_m0": 0.308,
    "omega_lambda0": 0.692,
    "h": 0.678,
}

# Prada12 的 Bolshoi/MultiDark 拟合所用 WMAP5-like cosmology。
PRADA_COSMOLOGY = {
    "name": "WMAP5-like (Prada12)",
    "omega_m0": 0.27,
    "omega_lambda0": 0.73,
    "h": 0.70,
}

# 球对称 top-hat 坍缩的线性临界密度。Ludlow16 Appendix C 取 1.686。
DELTA_SC = 1.686


def scale_factor(z):
    """把红移 z 转换为尺度因子 a=1/(1+z)。

    这是本阶段的完整示范函数。np.asarray 使函数可以同时接受一个数和数组。
    """
    # np.asarray 会把输入统一转换成 NumPy 浮点数组：
    #   scale_factor(1.0)                 可以接受单个数；
    #   scale_factor([0.0, 1.0, 2.0])     也可以接受一列数。
    # 后面的同一条公式会自动应用到每一个数组元素，这叫“向量化”。
    z = np.asarray(z, dtype=float)

    # np.any(condition) 检查数组中是否至少有一个元素满足 condition。
    # 这里主动拒绝负红移，能让输入错误尽早暴露，而不是产生难以理解的结果。
    if np.any(z < 0.0):
        raise ValueError("本次复现只使用 z >= 0。")

    # 注意：1.0 而不是 1 只是为了直观强调我们希望得到浮点数结果。
    return 1.0 / (1.0 + z)


def mass_to_hinv_msun_value(mass_msun, h):
    """将以 M_sun 给出的物理质量转换为以 h^-1 M_sun 表示的数值。

    若 M_phys = X h^-1 M_sun，则 M_phys = X/h M_sun，所以 X=h*M_phys。
    例如 1e12 M_sun 在 h=0.678 时写作 6.78e11 h^-1 M_sun。
    """
    mass_msun = np.asarray(mass_msun, dtype=float)
    if np.any(mass_msun <= 0.0):
        raise ValueError("halo mass 必须为正。")
    if h <= 0.0:
        raise ValueError("h 必须为正。")

    # 这里返回的只是“以 h^-1 M_sun 为单位时的数值”，不是给数组附加物理单位。
    # 例如 h=0.678 时：
    #   1.0e12 M_sun = 6.78e11 h^-1 M_sun。
    return h * mass_msun


def omega_m_z(z, omega_m0, omega_lambda0):
    """平直 matter+Lambda 宇宙在红移 z 的 Omega_m(z)。"""
    z = np.asarray(z, dtype=float)
    one_plus_z = 1.0 + z

    # E(z)^2 = [H(z)/H0]^2。这里不必先开方求 E，再把它平方。
    matter_term = omega_m0 * one_plus_z**3
    e_squared = matter_term + omega_lambda0
    return matter_term / e_squared

def omega_lambda_z(z, omega_m0, omega_lambda0):
    """平直 matter+Lambda 宇宙在红移 z 的 Omega_Lambda(z)。"""
    z = np.asarray(z, dtype=float)
    one_plus_z = 1.0 + z
    e_squared = omega_m0 * one_plus_z**3 + omega_lambda0
    return omega_lambda0 / e_squared

def growth_suppression(omega_m, omega_lambda):
    """Lahav/Carroll 近似中的增长抑制函数 g(Omega_m, Omega_Lambda)。
        g = (5 Omega_m / 2) /
            [Omega_m^(4/7) - Omega_Lambda
             + (1 + Omega_m/2)(1 + Omega_Lambda/70)].

    """
    omega_m = np.asarray(omega_m, dtype=float)
    omega_lambda = np.asarray(omega_lambda, dtype=float)
    denominator = (
        omega_m ** (4.0 / 7.0)
        - omega_lambda
        + (1.0 + omega_m / 2.0) * (1.0 + omega_lambda / 70.0)
    )
    numerator = 5.0 * omega_m / 2.0
    return numerator / denominator


def growth_factor_lahav(z, omega_m0, omega_lambda0):
    """返回归一化为 D(0)=1 的线性增长因子。

    TODO 3b：利用前面的三个函数完成

        D(z) = [g(z)/g(0)] / (1+z).

  """
    z = np.asarray(z, dtype=float)
    current_omega_m = omega_m_z(z, omega_m0, omega_lambda0)
    current_omega_lambda = omega_lambda_z(z, omega_m0, omega_lambda0)
    g_z = growth_suppression(current_omega_m, current_omega_lambda)
    g_0 = growth_suppression(omega_m0, omega_lambda0)
    unnormalized_growth = g_z / (1.0 + z)
    return unnormalized_growth / g_0


def ludlow_xi(mass_msun):
    """计算 Ludlow16 Appendix C9 的无量纲质量变量 xi。

    原始公式是

        xi = [M / (1e10 h^-1 M_sun)]^(-1).

    函数说明字符串必须紧跟在 def 后面；放在 return 后面时，Python 不会把它
    识别为该函数的文档。
    """
    h = LUDLOW_COSMOLOGY["h"]
    mass_msun = np.asarray(mass_msun, dtype=float)
    mass_hinv_msun = mass_to_hinv_msun_value(mass_msun, h)
    mass_ratio = mass_hinv_msun / 1.0e10
    return mass_ratio**(-1.0)


def sigma_ludlow16(mass_msun, z):
    """Ludlow16 Appendix C8 的线性 rms 涨落 sigma(M,z)。

    公式是

        sigma(M,z) = D(z) * 22.26 * xi^0.292
                     / [1 + 1.53*xi^0.275 + 3.36*xi^0.198].

    其中 D(z) 使用第一阶段完成的 Planck growth_factor_lahav。
    """
    z = np.asarray(z, dtype=float)
    mass_msun = np.asarray(mass_msun, dtype=float)
    xi = ludlow_xi(mass_msun)
    omega_m0 = LUDLOW_COSMOLOGY["omega_m0"]
    omega_lambda0 = LUDLOW_COSMOLOGY["omega_lambda0"]
    growth = growth_factor_lahav(z, omega_m0, omega_lambda0)
    numerator = 22.26 * xi**0.292
    denominator = 1.0 + 1.53 * xi**0.275 + 3.36 * xi**0.198
    return growth * numerator / denominator


def peak_height_ludlow16(mass_msun, z):
    """由 nu=delta_sc/sigma 计算 Ludlow16 使用的峰高。"""
    sigma = sigma_ludlow16(mass_msun, z)
    return DELTA_SC / sigma


def ludlow_c0(z):
    """计算 Ludlow16 Appendix C2 的浓度归一化参数 c0(z)。

    在完整的平滑破幂律中，浓度写成

        c(nu,z) = c0(z) * (nu/nu0)^(-gamma1)
                  * [1 + (nu/nu0)^(1/beta)]
                    ^[-beta*(gamma2-gamma1)].

    现在只实现第一个组成部分

        c0(z) = 3.395 * (1+z)^(-0.215).

    如何理解 c0：
    - 它控制整条 c(nu) 曲线的总体纵向高度，所以称为“归一化参数”；
    - 它随红移升高而缓慢减小；
    - 它本身不是最终浓度；
    - 严格来说，nu=nu0 时的最终浓度也不恰好等于 c0，因为方括号中的
      平滑转折因子此时不等于 1。

    c0 只依赖红移 z，与质量 M 和峰高 nu 都无关。
    """
    z = np.asarray(z, dtype=float)
    return 3.395 * (1.0 + z) ** (-0.215)


def ludlow_transition_beta(z):
    """返回平滑转折宽度 beta(z)，对应 Ludlow16 Appendix C3。

    beta 越大，两段幂律之间的转折越宽。这里把函数名写成
    ``transition_beta``，以免和后面的质量吸积史参数 beta 混淆。
    """
    z = np.asarray(z, dtype=float)
    return 0.307 * (1.0 + z) ** 0.540


def ludlow_gamma1(z):
    """返回低峰高端的渐近斜率参数 gamma1(z)，Appendix C4。"""
    z = np.asarray(z, dtype=float)
    return 0.628 * (1.0 + z) ** (-0.047)


def ludlow_gamma2(z):
    """返回高峰高端的渐近斜率参数 gamma2(z)，Appendix C5。"""
    z = np.asarray(z, dtype=float)
    return 0.317 * (1.0 + z) ** (-0.893)


def ludlow_nu0_published(z):
    """严格计算 Ludlow16 Appendix C6 印刷出的 nu0(z)。

    这个函数保留原始公式，便于研究拟合本身。需要注意一个真实的外推问题：
    原文虽然声明适用至大约 z=9，但印刷出的多项式在 z 约为 7.8 时穿过零，
    此后 nu0 为负，无法代入非整数幂。PBH 论文却把曲线画到了 z=12，因而
    必然还使用了原文没有说明的延拓处理。
    """
    z = np.asarray(z, dtype=float)
    a = scale_factor(z)

    polynomial = (
        4.135
        - 0.564 * a ** (-1.0)
        - 0.210 * a ** (-2.0)
        + 0.0557 * a ** (-3.0)
        - 0.00348 * a ** (-4.0)
    )

    growth = growth_factor_lahav(
        z,
        LUDLOW_COSMOLOGY["omega_m0"],
        LUDLOW_COSMOLOGY["omega_lambda0"],
    )
    return polynomial / growth


def ludlow_nu0(z):
    """返回用于 PBH Fig.1--2 的 nu0(z)，并显式处理高红移外推。

    在 z<=6 时使用出版公式；z>6 时从 z=6 的值按 (1+z)^(-1) 作连续延拓。
    这样会避开多项式接近零的病态区，并重现 PBH 论文图中不同质量曲线在
    高红移处继续平滑下降、逐渐靠近的行为。

    延拓的指数 -1 是根据 PBH Fig.1--2 的高红移趋势选取的复现约定，不是
    Ludlow16 发表的额外公式。后续进行数字化误差比较时必须把它列入误差账本。

    若研究原始拟合，请直接调用 ``ludlow_nu0_published``。
    """
    z = np.asarray(z, dtype=float)
    pivot_redshift = 6.0
    published_value = ludlow_nu0_published(np.minimum(z, pivot_redshift))
    continuation = ludlow_nu0_published(pivot_redshift) * (
        (1.0 + z) / (1.0 + pivot_redshift)
    ) ** (-1.0)
    return np.where(z <= pivot_redshift, published_value, continuation)


def concentration_ludlow16(mass_msun, z):
    """计算 Ludlow16 Appendix C1 的中位浓度 c_200(M,z)。

    为了看清平滑破幂律的结构，先定义无量纲比值

        ratio = nu(M,z) / nu0(z).

    第一项 ``ratio**(-gamma1)`` 给出低峰高一侧的趋势；方括号项让
    斜率在 ratio 经过 1 附近时平滑改变，而不是突然折断。
    """
    mass_msun = np.asarray(mass_msun, dtype=float)
    z = np.asarray(z, dtype=float)

    nu = peak_height_ludlow16(mass_msun, z)
    c0 = ludlow_c0(z)
    beta = ludlow_transition_beta(z)
    gamma1 = ludlow_gamma1(z)
    gamma2 = ludlow_gamma2(z)
    nu0 = ludlow_nu0(z)

    ratio = nu / nu0
    low_nu_power_law = ratio ** (-gamma1)
    smooth_transition = (
        1.0 + ratio ** (1.0 / beta)
    ) ** (-beta * (gamma2 - gamma1))
    return c0 * low_nu_power_law * smooth_transition


# -----------------------------------------------------------------------------
# Prada12：原始模型使用自己的时间变量 x 和未归一化到 D(0)=1 的增长解。
# -----------------------------------------------------------------------------


def prada_time_variable(z):
    """计算 Prada12 Eq.13 的时间变量 x。

    x 把尺度因子与 Omega_Lambda/Omega_m 组合起来。在 WMAP5 参数下，
    z=0 时 x=1.3931，正是 Prada12 在 B0、B1 中选取的参考点。
    """
    a = scale_factor(z)
    ratio = (
        PRADA_COSMOLOGY["omega_lambda0"]
        / PRADA_COSMOLOGY["omega_m0"]
    )
    return ratio ** (1.0 / 3.0) * a


def growth_factor_prada12(z):
    """计算 Prada12 Eq.12 的增长函数 D(a)。

    这里忠实使用 Prada12 与其 sigma 拟合配套的积分表达式，而不使用前面的
    Lahav 近似。Eq.12 右侧先给出早期满足 D/a -> 1 的增长解；Prada12 的
    文字定义要求最终增长因子在 z=0 归一化为 1，所以还要除以今天的值。

    ``quad`` 一次积分一个 x。下面先把任意形状数组拉平成一列，逐个积分后再
    恢复原形状，所以函数同时支持单个红移和红移数组。
    """
    x = np.asarray(prada_time_variable(z), dtype=float)

    def integrate_to_x(x_value):
        integrand = lambda u: u ** 1.5 / (1.0 + u**3) ** 1.5
        integral, _ = quad(integrand, 0.0, float(x_value))
        return integral

    integrals = np.array(
        [integrate_to_x(x_value) for x_value in x.ravel()]
    ).reshape(x.shape)

    omega_ratio = (
        PRADA_COSMOLOGY["omega_m0"]
        / PRADA_COSMOLOGY["omega_lambda0"]
    )
    prefactor = 2.5 * omega_ratio ** (1.0 / 3.0)
    unnormalized = (
        prefactor * np.sqrt(1.0 + x**3) / x**1.5 * integrals
    )

    x_today = float(prada_time_variable(0.0))
    integral_today = integrate_to_x(x_today)
    value_today = (
        prefactor
        * np.sqrt(1.0 + x_today**3)
        / x_today**1.5
        * integral_today
    )
    return unnormalized / value_today


def sigma_prada12(mass_msun, z):
    """计算 Prada12 Eq.23 的 sigma(M,z)。

    y 与 Ludlow16 的 xi 作用类似，也是无量纲反质量变量，但参考质量改为
    1e12 h^-1 M_sun，并使用 Prada12 的 h=0.70。
    """
    mass_hinv_msun = mass_to_hinv_msun_value(
        mass_msun, PRADA_COSMOLOGY["h"]
    )
    y = (mass_hinv_msun / 1.0e12) ** (-1.0)
    growth = growth_factor_prada12(z)
    numerator = 16.9 * y**0.41
    denominator = 1.0 + 1.102 * y**0.20 + 6.22 * y**0.333
    return growth * numerator / denominator


def prada_c_min(x):
    """Prada12 Eq.19：给定时间变量 x 时 U 形曲线的最小浓度。"""
    x = np.asarray(x, dtype=float)
    c0, c1 = 3.681, 5.033
    alpha, x0 = 6.948, 0.424
    transition = np.arctan(alpha * (x - x0)) / np.pi + 0.5
    return c0 + (c1 - c0) * transition


def prada_inverse_sigma_min(x):
    """Prada12 Eq.20：U 形最低点所对应的 sigma^{-1}。"""
    x = np.asarray(x, dtype=float)
    inverse_sigma0, inverse_sigma1 = 1.047, 1.646
    beta, x1 = 7.386, 0.526
    transition = np.arctan(beta * (x - x1)) / np.pi + 0.5
    return inverse_sigma0 + (inverse_sigma1 - inverse_sigma0) * transition


def prada_b0(x):
    """Prada12 Eq.18 的纵向重标度 B0(x)。"""
    return prada_c_min(x) / prada_c_min(1.393)


def prada_b1(x):
    """Prada12 Eq.18 的横向重标度 B1(x)。"""
    return prada_inverse_sigma_min(x) / prada_inverse_sigma_min(1.393)


def prada_universal_concentration(sigma_prime):
    """Prada12 Eqs.16--17 的近普适 U 形浓度函数。"""
    sigma_prime = np.asarray(sigma_prime, dtype=float)
    return (
        2.881
        * ((sigma_prime / 1.257) ** 1.022 + 1.0)
        * np.exp(0.060 / sigma_prime**2)
    )


def concentration_prada12(mass_msun, z, cap_high_peak=False):
    """计算 Prada12 Eqs.14--23 的浓度。

    ``cap_high_peak=False`` 返回 Prada12 原始 U 形拟合。

    PBH 论文指出 Prada12 在高红移、大峰高一侧会重新上翘，并在后续率计算中
    把该分支限制在相应红移的最小浓度。因此 ``cap_high_peak=True`` 时，仅对
    sigma_prime 小于 U 形最低点的位置使用 c_min(x)，低质量一侧不受影响。
    Fig.2 的论文曲线持续下降并在高红移趋于平台，故复现图使用这个选项。
    """
    x = prada_time_variable(z)
    sigma = sigma_prada12(mass_msun, z)
    sigma_prime = prada_b1(x) * sigma
    raw_concentration = prada_b0(x) * prada_universal_concentration(sigma_prime)

    if not cap_high_peak:
        return raw_concentration

    # B1 的定义使所有红移的最低点都映射到 sigma'=1/1.393。
    sigma_prime_at_minimum = 1.0 / 1.393
    return np.where(
        sigma_prime < sigma_prime_at_minimum,
        prada_c_min(x),
        raw_concentration,
    )


# -----------------------------------------------------------------------------
# PBH 论文 Appendix C：用今天的质量和浓度建立平均质量吸积史。
# -----------------------------------------------------------------------------


def nfw_g(concentration):
    """NFW 质量积分中反复出现的 g(c)=ln(1+c)-c/(1+c)。"""
    concentration = np.asarray(concentration, dtype=float)
    return np.log1p(concentration) - concentration / (1.0 + concentration)


def concentration_model(mass_msun, z, model, cap_prada=True):
    """用统一入口选择 Ludlow16 或 Prada12 浓度模型。"""
    model_name = model.lower()
    if model_name == "ludlow16":
        return concentration_ludlow16(mass_msun, z)
    if model_name == "prada12":
        return concentration_prada12(
            mass_msun, z, cap_high_peak=cap_prada
        )
    raise ValueError("model 必须是 'ludlow16' 或 'prada12'。")


def mass_accretion_parameters(mass0_msun, model):
    """由今天的质量 M0 求 Appendix C 的 z_-2、alpha 和 beta。

    注意这里的 alpha、beta 是质量吸积史参数，不是 Ludlow16 浓度公式中的
    gamma1 和 transition_beta。A_cosmo=798 是 PBH 论文给出的常数。
    """
    if model.lower() == "ludlow16":
        cosmology = LUDLOW_COSMOLOGY
    elif model.lower() == "prada12":
        cosmology = PRADA_COSMOLOGY
    else:
        raise ValueError("model 必须是 'ludlow16' 或 'prada12'。")

    concentration0 = concentration_model(
        mass0_msun, 0.0, model, cap_prada=False
    )
    omega_m0 = cosmology["omega_m0"]
    omega_lambda0 = cosmology["omega_lambda0"]
    a_cosmo = 798.0

    cube = (
        200.0
        * concentration0**3
        * nfw_g(1.0)
        / (a_cosmo * omega_m0 * nfw_g(concentration0))
        - omega_lambda0 / omega_m0
    )
    z_minus2 = np.cbrt(cube) - 1.0

    beta_mah = -3.0 / (1.0 + z_minus2)
    alpha_mah = (
        np.log(nfw_g(1.0) / nfw_g(concentration0))
        - beta_mah * z_minus2
    ) / np.log1p(z_minus2)
    return z_minus2, alpha_mah, beta_mah


def mass_accretion_history(mass0_msun, z, model):
    """计算平均主晕质量 M(z)=M0(1+z)^alpha exp(beta*z)。"""
    z = np.asarray(z, dtype=float)
    _, alpha_mah, beta_mah = mass_accretion_parameters(mass0_msun, model)
    return mass0_msun * (1.0 + z) ** alpha_mah * np.exp(beta_mah * z)


# -----------------------------------------------------------------------------
# Fig.1--2 的数据与绘图。数据函数和画图函数分开，方便后续直接复用曲线。
# -----------------------------------------------------------------------------


FIGURE_MASSES_MSUN = np.array([1.0e3, 1.0e6, 1.0e9, 1.0e12])
FIGURE_LINESTYLES = ["-", "--", "-.", ":"]


def make_figure1(output_directory):
    """生成论文 Fig.1：今天为 1e12 M_sun 晕的质量和浓度演化。"""
    output_directory = Path(output_directory)
    z = np.geomspace(0.1, 12.0, 400)
    mass = mass_accretion_history(1.0e12, z, "ludlow16")
    concentration = concentration_ludlow16(mass, z)

    figure, concentration_axis = plt.subplots(figsize=(7.0, 4.8))
    mass_axis = concentration_axis.twinx()

    concentration_axis.plot(
        z, concentration, color="tab:blue", linestyle="--", label=r"$C(z)$"
    )
    mass_axis.plot(z, mass, color="tab:red", label=r"$M(z)$")

    concentration_axis.set_xscale("log")
    mass_axis.set_yscale("log")
    concentration_axis.set_xlim(0.1, 12.0)
    concentration_axis.set_xlabel(r"Redshift $z$")
    concentration_axis.set_ylabel(r"Concentration $C$", color="tab:blue")
    mass_axis.set_ylabel(r"$M(z)\;[M_\odot]$", color="tab:red")
    concentration_axis.tick_params(axis="y", colors="tab:blue")
    mass_axis.tick_params(axis="y", colors="tab:red")

    lines = concentration_axis.lines + mass_axis.lines
    labels = [line.get_label() for line in lines]
    concentration_axis.legend(lines, labels, loc="upper right")
    figure.tight_layout()
    figure.savefig(output_directory / "fig1_reproduction.png", dpi=220)
    plt.close(figure)

    data = np.column_stack((z, mass, concentration))
    np.savetxt(
        output_directory / "fig1_reproduction.csv",
        data,
        delimiter=",",
        header="z,mass_msun,concentration_ludlow16",
        comments="",
    )
    return z, mass, concentration


def figure2_model_data(z, model):
    """返回 Fig.2 某个浓度模型的四条质量轨迹和浓度轨迹。"""
    mass_tracks = []
    concentration_tracks = []

    for mass0_msun in FIGURE_MASSES_MSUN:
        mass = mass_accretion_history(mass0_msun, z, model)
        concentration = concentration_model(
            mass, z, model, cap_prada=(model.lower() == "prada12")
        )
        mass_tracks.append(mass)
        concentration_tracks.append(concentration)

    return np.asarray(mass_tracks), np.asarray(concentration_tracks)


def make_figure2(output_directory):
    """生成论文 Fig.2：Ludlow16 与 Prada12 的浓度演化比较。"""
    output_directory = Path(output_directory)
    z = np.linspace(0.0, 12.0, 401)
    ludlow_mass, ludlow_concentration = figure2_model_data(z, "ludlow16")
    prada_mass, prada_concentration = figure2_model_data(z, "prada12")

    figure, axes = plt.subplots(
        2, 1, figsize=(7.2, 8.0), sharex=True, constrained_layout=True
    )

    for index, mass0_msun in enumerate(FIGURE_MASSES_MSUN):
        label = rf"$M_0=10^{{{int(np.log10(mass0_msun))}}}\,M_\odot$"
        axes[0].plot(
            z,
            ludlow_concentration[index],
            linestyle=FIGURE_LINESTYLES[index],
            label=label,
        )
        axes[1].plot(
            z,
            prada_concentration[index],
            linestyle=FIGURE_LINESTYLES[index],
            label=label,
        )

    axes[0].set_title("Ludlow16")
    axes[1].set_title("Prada12 — high-peak upturn capped")
    axes[0].set_ylim(2.5, 26.0)
    axes[1].set_ylim(3.0, 32.0)
    for axis in axes:
        axis.set_xlim(0.0, 12.0)
        axis.set_ylabel(r"Concentration $C$")
        axis.legend(loc="upper right")
    axes[1].set_xlabel(r"Redshift $z$")

    figure.savefig(output_directory / "fig2_reproduction.png", dpi=220)
    plt.close(figure)

    ludlow_table = np.column_stack((z, ludlow_mass.T, ludlow_concentration.T))
    prada_table = np.column_stack((z, prada_mass.T, prada_concentration.T))
    mass_columns = [f"mass_track_M0_{mass:.0e}_msun" for mass in FIGURE_MASSES_MSUN]
    concentration_columns = [
        f"concentration_M0_{mass:.0e}" for mass in FIGURE_MASSES_MSUN
    ]
    header = ",".join(["z", *mass_columns, *concentration_columns])
    np.savetxt(
        output_directory / "fig2_ludlow16.csv",
        ludlow_table,
        delimiter=",",
        header=header,
        comments="",
    )
    np.savetxt(
        output_directory / "fig2_prada12.csv",
        prada_table,
        delimiter=",",
        header=header,
        comments="",
    )
    return z, ludlow_concentration, prada_concentration


def print_key_values():
    """打印少量用于人工对照论文图像的关键数值，不执行自动判定。"""
    print("\nKey values for manual comparison")
    print("--------------------------------")
    print(f"Ludlow nu0(z=0)       = {float(ludlow_nu0(0.0)):.6f}")
    print(f"Prada x(z=0)          = {float(prada_time_variable(0.0)):.6f}")
    print(f"Prada D(z=0), normalized = {float(growth_factor_prada12(0.0)):.6f}")

    for model in ("ludlow16", "prada12"):
        concentrations = concentration_model(
            FIGURE_MASSES_MSUN,
            0.0,
            model,
            cap_prada=(model == "prada12"),
        )
        print(f"{model:9s} C(M0,z=0)  = {concentrations}")


def main():
    """生成今天阶段的两张图、三份曲线数据和人工对照数值。"""
    output_directory = Path(__file__).resolve().parent
    make_figure1(output_directory)
    make_figure2(output_directory)
    print_key_values()
    print("\nCreated fig1_reproduction.png and fig2_reproduction.png")
    print("Created fig1_reproduction.csv, fig2_ludlow16.csv, fig2_prada12.csv")

if __name__ == "__main__":
    main()
