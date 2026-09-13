"""NFW 暗物质晕结构：9 月 10 日阶段的逐步实现。

本模块接收前两阶段已经得到的晕质量、浓度和红移，逐步构造：

    (M, C, z) -> (R_vir, R_s, rho_s) -> rho_NFW(r) -> M(<r)

当前第一小步只实现论文 Eq. (2) 的 NFW 特征密度对比。密度剖面、内部质量、
PBH 速度分布和两体捕获率将在理解相应公式后逐个加入。

这里没有自动自检、断言或独立测试模块。
"""

import numpy as np

import haloconcentration as hc


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


# TODO（下一小步）：实现 rho_crit(z) 和 rho_s(C,z)。
# TODO（之后）：实现 rho_NFW(r)、M(<r)、特征速度和截断 Maxwell 分布。
