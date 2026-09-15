"""暗物质晕中的 PBH 两体引力波捕获。

本模块实现目标论文 Eq. (11) 的捕获截面、速度平均 ``<Sigma v>``，以及
Appendix A14 对任意 PBH 质量分布的每晕捕获率。单色和连续质量谱使用同一个
每晕率函数，区别只在传入的质量积分节点和权重。

速度分布沿用 ``halostructure.py`` 中归一化的三维各向同性密度 ``f_3(v)``：

    4*pi*integral v^2*f_3(v) dv = 1.

因此速度平均显式使用球壳测度 ``4*pi*v^2 dv``，不会重复或遗漏 Jacobian。
这里没有自动自检、断言或独立测试模块。
"""

from pathlib import Path

import matplotlib
import numpy as np
from scipy.integrate import quad

import haloconcentration as hc
import halostructure as hs
import pbhmassfunction as pmf


# 本脚本只把图片写入文件，不需要弹出交互式绘图窗口。显式选择 Agg 后端，
# 在 VS Code 终端或没有图形界面的环境中运行时也能稳定生成 PNG。
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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


def calculate_figure3_curves(halo_masses_msun=None):
    """计算论文 Fig. 3 的四条每晕直接捕获率曲线。

    Fig. 3 在红移 ``z=0`` 比较两个彼此独立的选择：

    1. 晕浓度模型：Ludlow16 或 Prada12；
    2. PBH 质量函数：30 M_sun 单色谱或 sigma=0.6 的对数正态谱。

    两个选择两两组合，恰好给出四条曲线。这里没有为单色谱另写捕获率
    公式；两种质量函数最终都被转换成 ``(质量节点, psi(m)dm 权重)``，然后
    传给同一个 ``capture_rate_per_halo_a14_per_year`` 函数。

    Parameters
    ----------
    halo_masses_msun : array-like, optional
        要计算的晕质量网格，单位为 M_sun。若省略，就在 10^3--10^15 M_sun
        之间使用对数均匀网格。对数网格意味着横轴每个十倍质量区间取得
        相同数量的计算点，适合论文中的双对数坐标。

    Returns
    -------
    halo_masses_msun : numpy.ndarray
        实际使用的一维晕质量网格。
    curves : dict[str, numpy.ndarray]
        四条曲线的数据。字典键同时标明浓度模型与 PBH 质量函数。
    """
    if halo_masses_msun is None:
        # 121 个点覆盖 12 个十倍区间，即每个 decade 有 10 个间隔。
        # 这已经足以画出平滑曲线，也不会让 Appendix A14 的积分过慢。
        halo_masses_msun = np.logspace(3.0, 15.0, 121)
    else:
        halo_masses_msun = np.asarray(halo_masses_msun, dtype=float)

    if halo_masses_msun.ndim != 1 or np.any(halo_masses_msun <= 0.0):
        raise ValueError("halo_masses_msun 必须是一维正数数组。")

    # 单色谱是一个离散质量节点：m=30 M_sun，质量分数权重为 1。
    monochromatic_nodes, monochromatic_weights = (
        pmf.monochromatic_mass_distribution(mass_msun=30.0)
    )

    # Fig. 3 使用 mu=ln(30 M_sun)、sigma=0.6 的对数正态谱。
    # pbhmassfunction.py 里的参数名 median_mass_msun 对应 exp(mu)=30 M_sun。
    def figure3_lognormal_pdf(mass_msun):
        return pmf.lognormal_mass_pdf(
            mass_msun,
            median_mass_msun=30.0,
            sigma=0.6,
        )

    # Gauss-Legendre 权重已经包含 psi(m_i)dm，因此后续双质量积分只需加权求和。
    lognormal_nodes, lognormal_weights = pmf.continuous_mass_quadrature(
        figure3_lognormal_pdf,
        point_count=96,
    )

    curves = {}
    for concentration_model in ("ludlow16", "prada12"):
        # capture_rate_per_halo_a14_per_year 当前接收单个晕质量，因此这里沿着
        # 晕质量网格逐点计算。每一点内部的 PBH 双质量积分仍由 NumPy 向量化。
        curves[f"{concentration_model}_monochromatic"] = np.array(
            [
                capture_rate_per_halo_a14_per_year(
                    halo_mass_msun,
                    z=0.0,
                    mass_nodes_msun=monochromatic_nodes,
                    mass_probability_weights=monochromatic_weights,
                    concentration_model=concentration_model,
                    f_pbh=1.0,
                )
                for halo_mass_msun in halo_masses_msun
            ]
        )
        curves[f"{concentration_model}_lognormal"] = np.array(
            [
                capture_rate_per_halo_a14_per_year(
                    halo_mass_msun,
                    z=0.0,
                    mass_nodes_msun=lognormal_nodes,
                    mass_probability_weights=lognormal_weights,
                    concentration_model=concentration_model,
                    f_pbh=1.0,
                )
                for halo_mass_msun in halo_masses_msun
            ]
        )

    return halo_masses_msun, curves


def make_figure3(output_directory="."):
    """绘制并保存 Fig. 3，同时导出构图所用的四列捕获率数据。

    生成文件：

    - ``fig3_reproduction.png``：四条曲线与低质量端局部放大图；
    - ``fig3_reproduction.csv``：晕质量和四条曲线的数值。

    返回值是图片路径与 CSV 路径，便于其他脚本以后调用本函数。
    """
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    halo_masses_msun, curves = calculate_figure3_curves()

    # 将曲线样式集中写在一个列表中。主图和插图循环使用同一列表，避免
    # 两处颜色、线型或标签不一致。
    curve_styles = [
        (
            "ludlow16_monochromatic",
            "Ludlow16 + mono.",
            "tab:blue",
            "-",
        ),
        (
            "prada12_monochromatic",
            "Prada12 + mono.",
            "tab:orange",
            "--",
        ),
        (
            "ludlow16_lognormal",
            "Ludlow16 + log-normal",
            "tab:green",
            "-.",
        ),
        (
            "prada12_lognormal",
            "Prada12 + log-normal",
            "tab:red",
            ":",
        ),
    ]

    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    for curve_key, label, color, line_style in curve_styles:
        axis.plot(
            halo_masses_msun,
            curves[curve_key],
            color=color,
            linestyle=line_style,
            linewidth=2.0,
            label=label,
        )

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(1.0e3, 1.0e15)
    axis.set_ylim(1.0e-15, 1.0e-10)
    axis.set_xlabel(r"Halo Mass ($M_\odot$)")
    axis.set_ylabel(r"$R_{\rm halo}\;(\mathrm{yr}^{-1})$")
    axis.legend(loc="upper left", fontsize=8.5, frameon=True)
    axis.grid(which="major", color="0.82", linewidth=0.7)
    axis.grid(which="minor", color="0.92", linewidth=0.45)

    # 论文原图在右下角放大低质量端。插图仍然使用双对数坐标，只缩小
    # 显示范围；它不是额外计算的一组数据。
    inset_axis = axis.inset_axes([0.57, 0.12, 0.40, 0.35])
    for curve_key, _, color, line_style in curve_styles:
        inset_axis.plot(
            halo_masses_msun,
            curves[curve_key],
            color=color,
            linestyle=line_style,
            linewidth=1.4,
        )
    inset_axis.set_xscale("log")
    inset_axis.set_yscale("log")
    inset_axis.set_xlim(1.0e3, 1.0e6)
    inset_axis.set_ylim(1.0e-15, 2.0e-13)
    inset_axis.tick_params(axis="both", which="both", labelsize=7)
    inset_axis.grid(which="major", color="0.84", linewidth=0.55)
    inset_axis.grid(which="minor", color="0.93", linewidth=0.35)

    figure.tight_layout()
    figure_path = output_directory / "fig3_reproduction.png"
    figure.savefig(figure_path, dpi=240, bbox_inches="tight")
    plt.close(figure)

    # CSV 的列顺序与图例顺序一致。第一行是列名，便于随后用 NumPy、
    # pandas 或电子表格软件读取并检查任意质量点。
    csv_path = output_directory / "fig3_reproduction.csv"
    table = np.column_stack(
        [
            halo_masses_msun,
            curves["ludlow16_monochromatic"],
            curves["prada12_monochromatic"],
            curves["ludlow16_lognormal"],
            curves["prada12_lognormal"],
        ]
    )
    np.savetxt(
        csv_path,
        table,
        delimiter=",",
        header=(
            "halo_mass_msun,"
            "ludlow16_monochromatic_per_year,"
            "prada12_monochromatic_per_year,"
            "ludlow16_lognormal_per_year,"
            "prada12_lognormal_per_year"
        ),
        comments="",
        fmt="%.10e",
    )

    return figure_path, csv_path


def main():
    """直接运行本脚本时生成 Fig. 3；被导入时不会自动绘图。"""
    figure_path, csv_path = make_figure3()
    print(f"Fig. 3 已保存到：{figure_path.resolve()}")
    print(f"Fig. 3 数据已保存到：{csv_path.resolve()}")


if __name__ == "__main__":
    main()
