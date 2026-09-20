"""不经过随机采样，直接由概率分布计算半长轴 ``a`` 的期望值。

这个脚本分别计算当前项目中的两种分布：

1. Appendix D 严格联合分布在
   ``10^-6 pc < a < 1 pc``、``0 < j < 1`` 上整体截断后的期望；
2. ``sample_figure19_orbital_parameters`` 使用的论文顺序有效分布的期望。

第一种情况先把联合密度对 ``j`` 积分成未归一化边缘密度
``P_tilde(a)``，然后计算

``E[a] = integral(a * P_tilde(a) da) / integral(P_tilde(a) da)``。

这里在 ``t = ln(a)`` 上做数值积分。这只是改善跨越六个数量级时的数值
分辨率；积分测度仍然是 ``da = exp(t) dt``，没有把物理分布改成
``d ln(a)``。

第二种情况按照论文顺序分布的定义，对 ``j`` 和
``X = (a / a_scale)^(3/4)`` 做确定性的 Gauss-Legendre 积分。它使用的有效
尺度来自 Fig. 19 经验复现分支，不能视为 Appendix D 唯一推出的理论参数。
"""

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

from p_a_j import (
    AU_PER_PC,
    DEFAULT_ALPHA,
    DEFAULT_RHO_EQ_MSUN_PC3,
    DEFAULT_SIGMA_EQ,
    FIGURE19_EFFECTIVE_ANGULAR_MOMENTUM_SCALE,
    FIGURE19_EFFECTIVE_SEMIMAJOR_AXIS_SCALE_PC,
    angular_momentum_cdf,
    semi_major_axis_marginal_pdf,
)


def expected_a_strict_appendix_d_pc(
    pbh_mass_msun=30.0,
    f_pbh=1.0,
    rho_eq_msun_pc3=DEFAULT_RHO_EQ_MSUN_PC3,
    sigma_eq=DEFAULT_SIGMA_EQ,
    alpha=DEFAULT_ALPHA,
    minimum_semi_major_axis_pc=1.0e-6,
    maximum_semi_major_axis_pc=1.0,
):
    """返回严格 Appendix D 截断联合分布的 ``E[a]``，单位为 pc。

    ``semi_major_axis_marginal_pdf`` 已经把联合密度在 ``0<j<1`` 上积分，
    但没有在有限 ``a`` 区间上归一化，所以期望必须除以归一化积分。
    """
    minimum_log_a = np.log(float(minimum_semi_major_axis_pc))
    maximum_log_a = np.log(float(maximum_semi_major_axis_pc))

    def marginal_density_at_log_a(log_a):
        semi_major_axis_pc = np.exp(log_a)
        return float(
            semi_major_axis_marginal_pdf(
                semi_major_axis_pc,
                pbh_mass_msun=pbh_mass_msun,
                f_pbh=f_pbh,
                rho_eq_msun_pc3=rho_eq_msun_pc3,
                sigma_eq=sigma_eq,
                alpha=alpha,
            )
        )

    # da = a d(ln a)，所以归一化积分中的 Jacobian 是 a。
    normalization = quad(
        lambda log_a: (
            marginal_density_at_log_a(log_a) * np.exp(log_a)
        ),
        minimum_log_a,
        maximum_log_a,
        epsabs=1.0e-12,
        epsrel=1.0e-10,
        limit=300,
    )[0]

    # 分子原本是 integral a*P(a) da；换元后多一个 da 的 a，故为 a^2。
    first_moment = quad(
        lambda log_a: (
            marginal_density_at_log_a(log_a) * np.exp(2.0 * log_a)
        ),
        minimum_log_a,
        maximum_log_a,
        epsabs=1.0e-12,
        epsrel=1.0e-10,
        limit=300,
    )[0]

    return first_moment / normalization


def expected_a_figure19_paper_order_pc(
    sigma_eq=DEFAULT_SIGMA_EQ,
    minimum_semi_major_axis_pc=1.0e-6,
    maximum_semi_major_axis_pc=1.0,
    angular_momentum_scale_for_sampling=(
        FIGURE19_EFFECTIVE_ANGULAR_MOMENTUM_SCALE
    ),
    semi_major_axis_scale_pc=(
        FIGURE19_EFFECTIVE_SEMIMAJOR_AXIS_SCALE_PC
    ),
    angular_momentum_node_count=128,
    transformed_axis_node_count=256,
):
    """返回论文顺序有效分布的 ``E[a]``，单位为 pc。

    对每个 ``j``，条件密度满足

    ``P(X|j) proportional to exp(-X) * g[j/(cX)]``，

    其中 ``g(y)=y^2/(1+y^2)^(3/2)``、
    ``c=0.5*sqrt(1+sigma_eq^2)``。程序先计算 ``E[a|j]``，再对截断后的
    固定尺度 ``P(j)`` 积分。整个过程不生成随机数。
    """
    j_nodes_raw, j_weights_raw = leggauss(
        int(angular_momentum_node_count)
    )
    j = 0.5 * (j_nodes_raw + 1.0)
    j_weights = 0.5 * j_weights_raw

    minimum_x = (
        minimum_semi_major_axis_pc / semi_major_axis_scale_pc
    ) ** 0.75
    maximum_x = (
        maximum_semi_major_axis_pc / semi_major_axis_scale_pc
    ) ** 0.75

    x_nodes_raw, x_weights_raw = leggauss(
        int(transformed_axis_node_count)
    )
    x = (
        0.5 * (maximum_x - minimum_x) * x_nodes_raw
        + 0.5 * (maximum_x + minimum_x)
    )
    x_weights = 0.5 * (maximum_x - minimum_x) * x_weights_raw

    # 论文顺序第一步使用固定有效尺度的 P(j)，并截断、归一化到 0<j<1。
    y_prior = j / angular_momentum_scale_for_sampling
    j_density = (
        y_prior**2
        / (j * (1.0 + y_prior**2) ** 1.5)
        / angular_momentum_cdf(
            1.0,
            angular_momentum_scale_for_sampling,
        )
    )

    # 对每一个 j 同时计算条件分布的归一化积分与 a 的一阶矩。
    scale_coefficient = 0.5 * np.sqrt(1.0 + sigma_eq**2)
    conditional_y = j[:, np.newaxis] / (
        scale_coefficient * x[np.newaxis, :]
    )
    conditional_shape = (
        conditional_y**2 / (1.0 + conditional_y**2) ** 1.5
    )
    unnormalized_x_density = (
        np.exp(-x)[np.newaxis, :] * conditional_shape
    )

    conditional_normalization = np.sum(
        unnormalized_x_density * x_weights,
        axis=1,
    )
    semi_major_axis_pc = (
        semi_major_axis_scale_pc * x ** (4.0 / 3.0)
    )
    conditional_first_moment = np.sum(
        unnormalized_x_density
        * semi_major_axis_pc[np.newaxis, :]
        * x_weights,
        axis=1,
    )
    conditional_expected_a = (
        conditional_first_moment / conditional_normalization
    )

    return np.sum(j_weights * j_density * conditional_expected_a)


if __name__ == "__main__":
    strict_mean_pc = expected_a_strict_appendix_d_pc()
    paper_order_mean_pc = expected_a_figure19_paper_order_pc()

    print("由概率密度直接积分得到的半长轴期望：")
    print(
        "严格 Appendix D 联合分布："
        f" {strict_mean_pc:.8f} pc"
        f" = {strict_mean_pc * AU_PER_PC:.2f} AU"
    )
    print(
        "论文顺序的 Fig. 19 有效分布："
        f" {paper_order_mean_pc:.8f} pc"
        f" = {paper_order_mean_pc * AU_PER_PC:.2f} AU"
    )
