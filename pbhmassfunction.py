"""PBH 质量分布：9 月 11 日阶段的逐步实现。

本模块只描述“PBH 总质量怎样分配到不同质量区间”，不负责 NFW 晕结构、
速度平均或两体捕获率。后续捕获率代码只通过这里公开的质量分布接口获取
``psi(m) dm``，从而可以替换质量分布而不修改晕模型。

当前实现论文 Appendix B2 的截断对数正态分布、Appendix B3 的 broken power-law、
Appendix B4 的 critical-collapse 分布，以及可由同一个 Appendix A14 积分器调用的
单色分布。Fig. 18 的具体参数、绘图和数据导出放在
``notebooks/fig18_pbh_mass_functions.ipynb``。

这里没有自动自检、断言或独立测试模块。
"""

import numpy as np
from scipy.special import gammainc, gammaln, ndtr


# 论文 Appendix B1 与 Fig. 18 指定的 PBH 质量积分范围，单位为 M_sun。
PBH_MASS_MIN_MSUN = 5.0
PBH_MASS_MAX_MSUN = 150.0


def lognormal_mass_pdf(
    mass_msun,
    median_mass_msun=30.0,
    sigma=0.71,
    mass_min_msun=PBH_MASS_MIN_MSUN,
    mass_max_msun=PBH_MASS_MAX_MSUN,
):
    """返回有限质量区间内归一化的对数正态 PBH 质量概率密度。

    论文 Appendix B2 的未截断表达式是

        psi_raw(m) = exp[-(ln(m)-mu)^2/(2*sigma^2)]
                     / (sqrt(2*pi)*sigma*m),

    其中 ``mu=ln(M_c)``，所以 ``M_c`` 是质量分布的中位数。论文随后要求
    所有分布在 5--150 M_sun 上归一化，因此这里还要除以该区间内的累计概率。

    Parameters
    ----------
    mass_msun : float 或 array-like
        要计算概率密度的 PBH 质量，单位为 M_sun。
    median_mass_msun : float, optional
        对数正态分布的中位质量 M_c，Fig. 18 默认使用 30 M_sun。
    sigma : float, optional
        ``ln(m)`` 的标准差，是无量纲的对数宽度。Fig. 18 默认使用 0.71。
    mass_min_msun, mass_max_msun : float, optional
        归一化和积分使用的质量下限、上限，默认分别为 5 和 150 M_sun。

    Returns
    -------
    psi : numpy.ndarray
        PBH 质量分数密度，单位为 M_sun^-1。区间外严格返回 0；在指定区间
        上对质量积分等于 1。

    Notes
    -----
    为了与 Appendix A2--A4 中 ``n_m=rho*f_PBH*psi(m)/m`` 一致，这里把
    ``psi(m) dm`` 解释为 PBH 总质量落入该质量区间的分数。它是关于线性质量
    测度 ``dm`` 的密度，所以带有 M_sun^-1 的单位；不要把它与关于
    ``d ln(m)`` 的密度混用。
    """
    mass_msun = np.asarray(mass_msun, dtype=float)

    if median_mass_msun <= 0.0:
        raise ValueError("median_mass_msun 必须为正。")
    if sigma <= 0.0:
        raise ValueError("sigma 必须为正。")
    if mass_min_msun <= 0.0 or mass_max_msun <= mass_min_msun:
        raise ValueError("质量上下限必须满足 0 < mass_min < mass_max。")

    mu = np.log(median_mass_msun)

    # ndtr(x) 是标准正态分布的累计分布函数 Phi(x)。变量替换
    # y=(ln(m)-mu)/sigma 后，有限区间内的原始概率就是 Phi(y_max)-Phi(y_min)。
    standardized_minimum = (np.log(mass_min_msun) - mu) / sigma
    standardized_maximum = (np.log(mass_max_msun) - mu) / sigma
    interval_probability = (
        ndtr(standardized_maximum) - ndtr(standardized_minimum)
    )

    inside_interval = (
        (mass_msun >= mass_min_msun) & (mass_msun <= mass_max_msun)
    )

    # 对区间外或非正质量，先用 1 M_sun 代替再计算对数，避免产生 log(0)
    # 或 log(负数) 的运行警告；最后仍由 np.where 将这些位置严格置零。
    safe_mass_msun = np.where(inside_interval, mass_msun, 1.0)
    log_distance = np.log(safe_mass_msun) - mu
    untruncated_density = np.exp(
        -(log_distance**2) / (2.0 * sigma**2)
    ) / (np.sqrt(2.0 * np.pi) * sigma * safe_mass_msun)

    return np.where(
        inside_interval,
        untruncated_density / interval_probability,
        0.0,
    )


def monochromatic_mass_distribution(mass_msun=30.0):
    """返回单色 PBH 质量分布的一点积分表示。

    连续记号中的单色质量函数是 Dirac 测度 ``delta(m-m0) dm``，不能在
    NumPy 网格上表示成一个普通的有限高度 PDF。对积分而言，它严格满足

        integral F(m) delta(m-m0) dm = F(m0).

    因此本函数返回一个质量节点 ``[m0]`` 和一个无量纲质量分数权重 ``[1]``。
    Appendix A14 的通用双质量积分器会照常对这些节点和权重求和，不需要另写
    一个“单色每晕捕获率”函数。
    """
    if mass_msun <= 0.0:
        raise ValueError("单色 PBH 质量必须为正。")

    return np.array([float(mass_msun)]), np.array([1.0])


def continuous_mass_quadrature(
    mass_pdf,
    mass_min_msun=PBH_MASS_MIN_MSUN,
    mass_max_msun=PBH_MASS_MAX_MSUN,
    point_count=96,
):
    """把连续质量密度 ``psi(m)`` 转成 Appendix A14 使用的节点和权重。

    这里采用 Gauss-Legendre 积分。返回的 ``probability_weights[i]`` 近似表示
    ``psi(m_i) dm``，所以双质量积分可写成两重加权求和。函数不会暗中重新
    归一化权重；质量 PDF 本身的归一化仍由各分布函数负责。
    """
    if mass_min_msun <= 0.0 or mass_max_msun <= mass_min_msun:
        raise ValueError("质量上下限必须满足 0 < mass_min < mass_max。")
    if point_count < 2:
        raise ValueError("point_count 必须至少为 2。")

    standard_nodes, standard_weights = np.polynomial.legendre.leggauss(
        point_count
    )
    half_width = 0.5 * (mass_max_msun - mass_min_msun)
    midpoint = 0.5 * (mass_max_msun + mass_min_msun)
    mass_nodes_msun = midpoint + half_width * standard_nodes
    probability_weights = (
        half_width
        * standard_weights
        * np.asarray(mass_pdf(mass_nodes_msun), dtype=float)
    )
    return mass_nodes_msun, probability_weights


def broken_power_law_mass_pdf(
    mass_msun,
    peak_mass_msun=31.1,
    alpha1=0.54,
    alpha2=5.6,
    mass_min_msun=PBH_MASS_MIN_MSUN,
    mass_max_msun=PBH_MASS_MAX_MSUN,
):
    """返回有限质量区间内归一化的 broken power-law PBH 质量密度。

    Appendix B3 在峰值质量 ``m_*`` 两侧分别取
    ``(m/m_*)**alpha1`` 和 ``(m/m_*)**(-alpha2)``。论文的 Fig. 18
    使用 ``m_*=31.1 M_sun``、``alpha1=0.54``、``alpha2=5.6``，并将分布
    在 5--150 M_sun 上重新归一化。

    Parameters
    ----------
    mass_msun : float 或 array-like
        要计算概率密度的 PBH 质量，单位为 M_sun。
    peak_mass_msun : float, optional
        分段幂律的峰值质量 ``m_*``，单位为 M_sun。
    alpha1, alpha2 : float, optional
        低质量端和高质量端的幂律指数；Appendix B3 要求 ``alpha1 > 0``、
        ``alpha2 > 1``。
    mass_min_msun, mass_max_msun : float, optional
        归一化和积分使用的质量下限、上限。

    Returns
    -------
    psi : numpy.ndarray
        关于线性质量测度 ``dm`` 的质量分数密度，单位为 M_sun^-1；区间外
        严格返回 0。
    """
    mass_msun = np.asarray(mass_msun, dtype=float)

    if peak_mass_msun <= 0.0:
        raise ValueError("peak_mass_msun 必须为正。")
    if alpha1 <= 0.0:
        raise ValueError("alpha1 必须为正。")
    if alpha2 <= 1.0:
        raise ValueError("alpha2 必须大于 1。")
    if mass_min_msun <= 0.0 or mass_max_msun <= mass_min_msun:
        raise ValueError("质量上下限必须满足 0 < mass_min < mass_max。")

    normalization = 0.0
    if mass_min_msun < peak_mass_msun:
        lower_upper = min(mass_max_msun, peak_mass_msun)
        normalization += peak_mass_msun * (
            (lower_upper / peak_mass_msun) ** (alpha1 + 1.0)
            - (mass_min_msun / peak_mass_msun) ** (alpha1 + 1.0)
        ) / (alpha1 + 1.0)
    if mass_max_msun > peak_mass_msun:
        upper_lower = max(mass_min_msun, peak_mass_msun)
        normalization += peak_mass_msun * (
            (upper_lower / peak_mass_msun) ** (1.0 - alpha2)
            - (mass_max_msun / peak_mass_msun) ** (1.0 - alpha2)
        ) / (alpha2 - 1.0)

    inside_interval = (
        (mass_msun >= mass_min_msun) & (mass_msun <= mass_max_msun)
    )
    safe_mass_ratio = np.where(
        inside_interval, mass_msun / peak_mass_msun, 1.0
    )
    unnormalized_density = np.where(
        safe_mass_ratio <= 1.0,
        safe_mass_ratio**alpha1,
        safe_mass_ratio ** (-alpha2),
    )
    return np.where(
        inside_interval, unnormalized_density / normalization, 0.0
    )


def critical_collapse_mass_pdf(
    mass_msun,
    mass_scale_msun=10.8,
    alpha=1.1,
    mass_min_msun=PBH_MASS_MIN_MSUN,
    mass_max_msun=PBH_MASS_MAX_MSUN,
):
    """返回有限质量区间内归一化的 critical-collapse PBH 质量密度。

    Appendix B4 定义的未截断分布为

        psi_raw(m) = alpha**2 * m**alpha
                     * exp[-(m/M_f)**alpha]
                     / (M_f**(1+alpha) * Gamma(1/alpha)).

    ``M_f`` 是临界坍缩的质量尺度。Fig. 18 图注使用
    ``M_f=10.8 M_sun``、``alpha=1.1``，而 Appendix B4 正文写 ``alpha=1``；
    这里采用图注值作为默认值，并将分布在论文指定的 5--150 M_sun 上归一化。

    Parameters
    ----------
    mass_msun : float 或 array-like
        要计算概率密度的 PBH 质量，单位为 M_sun。
    mass_scale_msun : float, optional
        临界坍缩质量尺度 ``M_f``，单位为 M_sun。
    alpha : float, optional
        临界坍缩指数，必须为正。
    mass_min_msun, mass_max_msun : float, optional
        归一化和积分使用的质量下限、上限。

    Returns
    -------
    psi : numpy.ndarray
        关于线性质量测度 ``dm`` 的质量分数密度，单位为 M_sun^-1；区间外
        严格返回 0。
    """
    mass_msun = np.asarray(mass_msun, dtype=float)

    if mass_scale_msun <= 0.0:
        raise ValueError("mass_scale_msun 必须为正。")
    if alpha <= 0.0:
        raise ValueError("alpha 必须为正。")
    if mass_min_msun <= 0.0 or mass_max_msun <= mass_min_msun:
        raise ValueError("质量上下限必须满足 0 < mass_min < mass_max。")

    gamma_shape = 1.0 + 1.0 / alpha
    transformed_minimum = (mass_min_msun / mass_scale_msun) ** alpha
    transformed_maximum = (mass_max_msun / mass_scale_msun) ** alpha
    interval_probability = (
        gammainc(gamma_shape, transformed_maximum)
        - gammainc(gamma_shape, transformed_minimum)
    )

    inside_interval = (
        (mass_msun >= mass_min_msun) & (mass_msun <= mass_max_msun)
    )
    safe_mass_ratio = np.where(
        inside_interval, mass_msun / mass_scale_msun, 1.0
    )
    transformed_mass = safe_mass_ratio**alpha
    log_untruncated_density = (
        2.0 * np.log(alpha)
        - np.log(mass_scale_msun)
        - gammaln(1.0 / alpha)
        + alpha * np.log(safe_mass_ratio)
        - transformed_mass
    )
    untruncated_density = np.exp(log_untruncated_density)

    return np.where(
        inside_interval,
        untruncated_density / interval_probability,
        0.0,
    )
