"""暗物质晕浓度模型与平均质量吸积史。

学习约定
--------
这个文件最初采用逐函数填空的学习方式。为了在 2026-09-20 前完成论文的
全流程精读与分层复现，从本版本开始改为“模块交付 + 逐段精讲”：函数保持
短小、注释解释公式和单位，同时每个阶段直接形成能够生成论文图像的程序。

当前版本提供：
1. Ludlow16 Appendix C 的浓度拟合；
2. Prada12 Eq. (23) 闭式近似，以及用 WMAP5 功率谱实现 Eqs. (12)--(22)
   的浓度计算；
3. PBH 论文 Appendix C 的平均质量吸积史。

Fig. 1 和 Fig. 2 的具体质量、红移网格、数据组织、绘图与导出分别放在
``notebooks/fig01_halo_mass_concentration_history.ipynb`` 和
``notebooks/fig02_concentration_model_comparison.ipynb``。

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

import numpy as np
from astropy.cosmology import WMAP5
from hmf import MassFunction
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

# 用 Prada12 Eqs. (12)--(22) 直接计算浓度时，需要由 WMAP5 线性功率谱得到
# sigma(M,z)。下面是 Bolshoi/MultiDark 所采用的归一化和原初谱指数。
PRADA_HMF_SIGMA_8 = 0.817
PRADA_HMF_SPECTRAL_INDEX = 0.962

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
    """
    z = np.asarray(z, dtype=float)
    return 3.395 * (1.0 + z) ** (-0.215)


def ludlow_transition_beta(z):
    """返回平滑转折宽度 beta(z)，对应 Ludlow16 Appendix C3。
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
    此后 nu0 为负，无法代入非整数幂。PBH 论文的浓度曲线画到了 z=12，
    但具体处理尚未确认，不能据此断言作者必然使用了某种延拓。
    该问题与 Fig. 9 的高质量吸积参数差异分开核对：后者只使用 C(M0,0)。
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


def linear_sigma_hmfcalc_wmap5(mass_msun, z, dlog10m=0.02):
    """用 ``hmf`` 的 WMAP5 线性功率谱计算 ``sigma(M,z)``。

    ``MassFunction.sigma`` 与所选的 halo mass-function 拟合式无关；这里借用
    同一个后端完成功率谱、top-hat 滤波和线性增长。公开质量输入为物理
    ``M_sun``，进入 ``hmf`` 前乘以 ``h``，转换成其 ``M_sun/h`` 数值。

    这个函数服务于 Prada12 Eqs. (12)--(22)。本文件中的
    :func:`sigma_prada12` 是 Eq. (23) 的闭式近似，两条路径同时保留。
    """
    mass_msun, z = np.broadcast_arrays(
        np.asarray(mass_msun, dtype=float),
        np.asarray(z, dtype=float),
    )
    dlog10m = float(dlog10m)

    if np.any(mass_msun <= 0.0):
        raise ValueError("计算 sigma 的 halo mass 必须为正。")
    if np.any(z < 0.0):
        raise ValueError("红移 z 不能为负。")
    if dlog10m <= 0.0:
        raise ValueError("dlog10m 必须为正。")

    hubble_h = float(WMAP5.h)
    minimum_mass_hinv = float(np.min(mass_msun)) * hubble_h
    maximum_mass_hinv = float(np.max(mass_msun)) * hubble_h

    # 在所需质量区间两端各多留 0.1 dex，使后面的对数插值不会落在边界外。
    sigma_calculator = MassFunction(
        Mmin=np.log10(minimum_mass_hinv) - 0.1,
        Mmax=np.log10(maximum_mass_hinv) + 0.1,
        dlog10m=dlog10m,
        z=0.0,
        hmf_model="PS",
        delta_c=DELTA_SC,
        cosmo_model=WMAP5,
        sigma_8=PRADA_HMF_SIGMA_8,
        n=PRADA_HMF_SPECTRAL_INDEX,
        transfer_model="CAMB",
        transfer_params={"extrapolate_with_eh": True},
    )

    physical_mass_grid = np.asarray(sigma_calculator.m, dtype=float) / hubble_h
    sigma0_grid = np.asarray(sigma_calculator.sigma, dtype=float)
    sigma0 = np.exp(
        np.interp(
            np.log(mass_msun),
            np.log(physical_mass_grid),
            np.log(sigma0_grid),
        )
    )
    # hmf 的增长因子接口在数组输入的 z=0 端点会经过样条插值，产生约 1e-4
    # 的端点误差。逐个标量调用只涉及已建好的增长样条，不会重新计算功率谱，
    # 同时与原先“每次传入一个红移”的结果保持一致。
    growth = np.asarray(
        [
            float(
                np.asarray(
                    sigma_calculator.growth.growth_factor(float(z_value)),
                    dtype=float,
                ).ravel()[0]
            )
            for z_value in z.ravel()
        ],
        dtype=float,
    ).reshape(z.shape)
    return sigma0 * growth


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

    # 1.393 是参考时间变量 x_ref，不是最低点的 sigma^{-1}。
    # B1 的定义使各红移的最低点映射到
    # sigma'_min = 1 / sigma_min^{-1}(x_ref)。
    sigma_prime_at_minimum = 1.0 / prada_inverse_sigma_min(1.393)
    return np.where(
        sigma_prime < sigma_prime_at_minimum,
        prada_c_min(x),
        raw_concentration,
    )


def concentration_prada12_hmf_sigma(
    mass_msun,
    z,
    cap_high_peak=True,
):
    """用 WMAP5 功率谱的 ``sigma(M,z)`` 计算 Prada12 浓度。

    组合公式仍是 Prada12 Eqs. (14)--(22)；与
    :func:`concentration_prada12` 的区别只在于这里不使用 Eq. (23) 的
    ``sigma`` 闭式近似。PBH 论文的率计算限制 Prada 高峰高上翘分支，所以
    ``cap_high_peak`` 默认开启。
    """
    mass_msun, z = np.broadcast_arrays(
        np.asarray(mass_msun, dtype=float),
        np.asarray(z, dtype=float),
    )
    sigma = linear_sigma_hmfcalc_wmap5(mass_msun, z)
    x = prada_time_variable(z)
    sigma_prime = prada_b1(x) * sigma
    raw_concentration = prada_b0(x) * prada_universal_concentration(
        sigma_prime
    )

    if not cap_high_peak:
        return raw_concentration

    sigma_prime_at_minimum = 1.0 / prada_inverse_sigma_min(1.393)
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
    """统一选择 Ludlow16、Prada12 闭式或 Prada12-HMF 浓度。"""
    model_name = model.lower()
    if model_name == "ludlow16":
        return concentration_ludlow16(mass_msun, z)
    if model_name == "prada12":
        return concentration_prada12(
            mass_msun, z, cap_high_peak=cap_prada
        )
    if model_name == "prada12_hmf_sigma":
        return concentration_prada12_hmf_sigma(
            mass_msun,
            z,
            cap_high_peak=cap_prada,
        )
    raise ValueError(
        "model 必须是 'ludlow16'、'prada12' 或 "
        "'prada12_hmf_sigma'。"
    )


def mass_accretion_parameters(mass0_msun, model):
    """由今天的质量 M0 求 Appendix C 的 z_-2、alpha 和 beta。

    注意这里的 alpha、beta 是质量吸积史参数，不是 Ludlow16 浓度公式中的
    gamma1 和 transition_beta。A_cosmo=798 是 PBH 论文给出的常数。
    """
    # Fig. 9 核对记录（2026-09-13）：原图 M0=1e13、1e14、1e15 M_sun 的
    # (alpha,beta) 分别为 (0.2702,-0.904)、(0.3312,-1.014)、(0.3922,-1.124)。
    # 它们来自 arXiv v1/v2 的 MAH_corea.pdf 矢量路径反求；低/高红移段结果一致。
    # 在本代码的 Planck 参数和 A=798 下，这些参数不能同时满足 Appendix C1-C3。
    # 推测：作者可能对高质量端另行指定/延拓了参数；尚无生成代码或作者说明，
    # 不能排除参数表、实现差异或图文不一致，也未发现高红移处切换公式的证据。
    # 本函数保留印刷公式；不将反求参数或跨轨迹排序作为修正。证据见
    # SIMULATING_PBH_MERGERS_精读与复现.md 第 7.4.1 节。
    if model.lower() == "ludlow16":
        cosmology = LUDLOW_COSMOLOGY
    elif model.lower() in {"prada12", "prada12_hmf_sigma"}:
        cosmology = PRADA_COSMOLOGY
    else:
        raise ValueError(
            "model 必须是 'ludlow16'、'prada12' 或 "
            "'prada12_hmf_sigma'。"
        )

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
