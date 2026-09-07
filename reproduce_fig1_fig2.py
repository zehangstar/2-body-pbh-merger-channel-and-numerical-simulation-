"""逐步复现 Aljaf & Cholis (2025) 的 Fig. 1 和 Fig. 2。

学习约定
--------
这个文件会分阶段扩展，而不是一次放入全部答案。

已完成阶段：1. 共同宇宙学基础；2A. Ludlow16 的 xi、sigma 与峰高 nu
当前阶段：2B-1. Ludlow16 平滑破幂律的归一化参数 c0(z)
下一阶段：2B-2. beta(z)、gamma1(z)、gamma2(z) 和 nu0(z)

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
5. 本文件中的 TODO 是练习位置。每次只完成一个 TODO，再运行脚本观察第一处
   报错是否移动到下一个 TODO。
"""

from __future__ import annotations

import numpy as np


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

    TODO 6a：请仿照 scale_factor 函数完成以下两步。

    第一步：用 np.asarray 把 z 转换成浮点数组；
    第二步：按照上式返回 c0(z)。Python 的乘方运算符是 **。

    暂时不要在这里调用 peak_height_ludlow16，因为 c0 只依赖红移 z，
    与质量 M 和峰高 nu 都无关。
    """
    # 在下一行写第一步，例如参考 scale_factor 中处理 z 的方式。

    # 在下一行写第二步。注意指数 -0.215 必须放在乘方运算的右边。

    raise NotImplementedError("TODO 6a：请实现 ludlow_c0(z)。")

def check_stage1():
    """完成 TODO 1--3 后运行的第一阶段自检。

    assert 的意思是“我确信后面的条件必须成立”。条件为 False 时程序立即停止，
    从而指出哪条物理或数学性质没有满足。全部 assert 都通过才会打印成功信息。
    """
    for cosmology in (LUDLOW_COSMOLOGY, PRADA_COSMOLOGY):
        # 字典用方括号和键名取值。例如 cosmology["omega_m0"] 会依次取得
        # 0.308 和 0.27，因为 for 循环会先后检查两套宇宙学。
        om0 = cosmology["omega_m0"]
        ol0 = cosmology["omega_lambda0"]

        # z=0 时应恢复输入的当前密度参数。
        assert np.isclose(omega_m_z(0.0, om0, ol0), om0)
        assert np.isclose(omega_lambda_z(0.0, om0, ol0), ol0)

        # np.array 建立一个包含四个测试红移的一维数组。
        # 平直宇宙中，Omega_m 与 Omega_Lambda 在任意 z 相加为 1。
        z_grid = np.array([0.0, 1.0, 3.0, 10.0])
        omega_sum = (
            omega_m_z(z_grid, om0, ol0)
            + omega_lambda_z(z_grid, om0, ol0)
        )
        assert np.allclose(omega_sum, 1.0)

        # 归一化和单调性：红移越高，线性增长因子越小。
        # np.diff(growth) 返回相邻元素之差；若严格递减，每个差都应小于 0。
        growth = growth_factor_lahav(z_grid, om0, ol0)
        assert np.isclose(growth[0], 1.0)
        assert np.all(np.isfinite(growth))
        assert np.all(growth > 0.0)
        assert np.all(np.diff(growth) < 0.0)

        print(f"{cosmology['name']}: D(z) = {growth}")

    converted_mass = mass_to_hinv_msun_value(1.0e12, LUDLOW_COSMOLOGY["h"])
    assert np.isclose(converted_mass, 6.78e11)
    print("Stage 1 checks passed.")


def check_stage2():
    """依次检查 Ludlow16 的 xi、sigma 和 nu。

    这些检查既验证数值，也验证物理趋势。完成一个 TODO 后再次运行脚本：
    如果报错移动到下一个 TODO，就说明刚完成的函数已经通过了前面的检查。
    """
    masses = np.array([1.0e3, 1.0e6, 1.0e9, 1.0e12])

    # C9：xi 与 M 成反比。第一个精确数值检查专门捕捉 h 单位转换错误。
    xi = ludlow_xi(masses)
    expected_xi_at_1e12 = 1.0 / 67.8
    assert np.isclose(xi[-1], expected_xi_at_1e12)
    assert np.all(np.isfinite(xi))
    assert np.all(xi > 0.0)
    assert np.all(np.diff(xi) < 0.0)
    print(f"Ludlow16 xi(M) = {xi}")

    # C8：固定 z=0 时，sigma 随质量增大而减小。
    sigma_z0 = sigma_ludlow16(masses, 0.0)
    assert np.isclose(sigma_z0[-1], 2.2119621585650333)
    assert np.all(np.isfinite(sigma_z0))
    assert np.all(sigma_z0 > 0.0)
    assert np.all(np.diff(sigma_z0) < 0.0)

    # 固定 M=1e12 M_sun 时，sigma 随红移升高而减小。
    z_grid = np.array([0.0, 1.0, 3.0, 10.0])
    sigma_redshift = sigma_ludlow16(1.0e12, z_grid)
    assert np.all(np.diff(sigma_redshift) < 0.0)
    print(f"Ludlow16 sigma(M, z=0) = {sigma_z0}")
    print(f"Ludlow16 sigma(1e12 M_sun, z) = {sigma_redshift}")

    # 峰高与 sigma 成反比，所以固定 z 时随质量增大，固定 M 时随红移增大。
    nu_z0 = peak_height_ludlow16(masses, 0.0)
    nu_redshift = peak_height_ludlow16(1.0e12, z_grid)
    assert np.allclose(nu_z0, DELTA_SC / sigma_z0)
    assert np.all(np.diff(nu_z0) > 0.0)
    assert np.all(np.diff(nu_redshift) > 0.0)
    print(f"Ludlow16 nu(M, z=0) = {nu_z0}")
    print(f"Ludlow16 nu(1e12 M_sun, z) = {nu_redshift}")
    print("Stage 2A checks passed.")


if __name__ == "__main__":
    check_stage1()
    check_stage2()
