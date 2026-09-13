"""暗物质晕质量吸积史：论文 Fig. 9 的逐步复现脚本。

这个文件对应 9 月 9 日的学习任务。它不重复实现 Ludlow16 浓度公式，
而是导入 ``haloconcentration.py`` 中已经完成并核对过的函数。

当前阶段完成三件事：

1. 为 Fig. 9 的一组现今晕质量计算 M(z)；
2. 使用该数据矩阵绘制并保存 Fig. 9；
3. 计算同一批晕的 C(z)、R_vir(z) 与 R_s(z)。

这里没有自动自检、断言或测试模块。
"""

from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import haloconcentration as hc


# Fig. 9 展示的“现今质量”范围是 10^3--10^15 M_sun。
# np.logspace(3, 15, 13) 会依次生成 10^3、10^4、...、10^15，
# 因而一共对应论文图例中的 13 条质量吸积轨迹。
FIG9_PRESENT_DAY_MASSES_MSUN = np.logspace(3.0, 15.0, 13)

# G 的单位写成 kpc (km/s)^2 M_sun^-1，使下面算出的半径直接以 kpc 表示。
GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN = 4.30091e-6
VIRIAL_OVERDENSITY = 200.0


def calculate_mass_accretion_tracks(
    z, present_day_masses_msun=FIG9_PRESENT_DAY_MASSES_MSUN
):
    """计算 Fig. 9 中不同现今质量晕的平均质量吸积轨迹。

    Parameters
    ----------
    z : float 或 array-like
        需要回溯的红移。z=0 表示今天；z 越大，表示回溯到越早的宇宙。
    present_day_masses_msun : array-like
        每条轨迹在 z=0 时的晕质量 M0，单位统一为太阳质量 M_sun。

    Returns
    -------
    mass_tracks_msun : numpy.ndarray
        二维数组，形状为 ``(现今质量的数量, 红移点的数量)``。
        第 i 行表示第 i 个 M0 对应的完整 M(z) 轨迹；
        第 j 列表示所有晕在同一个红移 z[j] 时的质量。

    Notes
    -----
    本文 Fig. 9 明确采用 Ludlow16。因此这里把 ``model`` 固定为
    ``"ludlow16"``，而不是让调用者在这一阶段任意切换模型。

    真正的质量吸积公式仍由 ``hc.mass_accretion_history`` 计算：

        M(z) = M0 * (1 + z)**alpha * exp(beta * z)

    其中每个 M0 都有自己的一组 z_-2、alpha 和 beta。
    """
    # 统一转成一维浮点数组。即使只输入一个红移或一个质量，
    # np.atleast_1d 也能让下面的循环和输出形状保持一致。
    z = np.atleast_1d(np.asarray(z, dtype=float))
    present_day_masses_msun = np.atleast_1d(
        np.asarray(present_day_masses_msun, dtype=float)
    )

    # 先批量求每个 M0 专属的 z_-2、alpha、beta，再用二维广播一次生成
    # 所有轨迹；行对应 M0，列对应 z。
    # 已知复现差异：当前公式的 1e15/1e14、1e14/1e13 M_sun 轨迹分别在
    # z≈7.464、10.795 相交，原 Fig. 9 对应曲线则不相交；不能视为已精确复现。
    # 原图参数提示高质量端可能有未说明的参数处理，但具体机制尚未确认。
    # 每一行必须保留同一 M0 的轨迹身份；逐红移排序会交换身份、掩盖差异。
    # 反求参数、公式一致性检查和推测边界见精读笔记第 7.4.1 节。
    _, alpha_mah, beta_mah = hc.mass_accretion_parameters(
        present_day_masses_msun, model="ludlow16"
    )
    mass0 = present_day_masses_msun[:, np.newaxis]
    return (
        mass0
        * (1.0 + z[np.newaxis, :]) ** alpha_mah[:, np.newaxis]
        * np.exp(beta_mah[:, np.newaxis] * z[np.newaxis, :])
    )


def make_figure9(output_directory):
    """使用质量吸积轨迹矩阵绘制并保存论文 Fig. 9 及对应数据。"""
    output_directory = Path(output_directory)

    # 原图的矢量路径包含 z=0 加上 49 个等间距的正红移点，范围到 z=12。
    # z=0 对确认 M(0)=M0 很重要，因此保留在 CSV 中；但 log(0) 不存在，
    # 所以绘图时只使用 positive_redshift 选出的 49 个正红移点。
    z = np.linspace(0.0, 12.0, 50)
    mass_tracks_msun = calculate_mass_accretion_tracks(z)
    positive_redshift = z > 0.0

    figure, axis = plt.subplots(figsize=(7.0, 4.8))
    line_styles = ("-", "--", "-.", ":")
    for index, mass0_msun in enumerate(FIG9_PRESENT_DAY_MASSES_MSUN):
        exponent = int(np.log10(mass0_msun))
        axis.plot(
            z[positive_redshift],
            mass_tracks_msun[index, positive_redshift],
            linestyle=line_styles[index % len(line_styles)],
            label=rf"$M=10^{{{exponent}}}M_\odot$",
        )

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(z[1], 12.0)
    axis.set_ylim(1.0e1, 2.0e15)
    axis.set_yticks(10.0 ** np.arange(3.0, 16.0, 3.0))
    axis.set_xlabel(r"Redshift $z$")
    axis.set_ylabel(r"$M(z)\;[M_\odot]$")
    axis.legend(
        loc="upper right",
        fontsize=6.5,
        handlelength=2.7,
        labelspacing=0.25,
    )
    figure.tight_layout()
    figure.savefig(output_directory / "fig9_reproduction.png", dpi=220)
    plt.close(figure)

    # CSV 的第一列是 z，后面 13 列依次对应 M0=10^3--10^15 M_sun。
    # 转置 mass_tracks_msun 后，每一行才会对应同一个红移，便于逐行读取。
    mass_columns = [
        f"mass_M0_{mass0_msun:.0e}_msun"
        for mass0_msun in FIG9_PRESENT_DAY_MASSES_MSUN
    ]
    np.savetxt(
        output_directory / "fig9_reproduction.csv",
        np.column_stack((z, mass_tracks_msun.T)),
        delimiter=",",
        header=",".join(["z", *mass_columns]),
        comments="",
    )
    return z, mass_tracks_msun


def calculate_halo_structure_tracks(
    z, present_day_masses_msun=FIG9_PRESENT_DAY_MASSES_MSUN
):
    """计算 ``M(z) -> C(z) -> R_vir(z) -> R_s(z)`` 的完整数据流。

    四个返回数组的形状均为 ``(现今质量的数量, 红移点的数量)``；质量单位
    为 M_sun，两个半径的单位为 kpc。virial 半径采用 M_200c 定义，即晕内
    平均密度为同一红移临界密度的 200 倍。
    """
    z = np.atleast_1d(np.asarray(z, dtype=float))
    mass_tracks_msun = calculate_mass_accretion_tracks(
        z, present_day_masses_msun
    )

    redshift_grid = z[np.newaxis, :]
    concentration_tracks = hc.concentration_model(
        mass_tracks_msun, redshift_grid, model="ludlow16"
    )

    cosmology = hc.LUDLOW_COSMOLOGY

    # H0=100 h km/s/Mpc。因为 1 Mpc=1000 kpc，所以换成 km/s/kpc 后
    # H0=0.1 h；再乘 E(z) 得到同一红移的 H(z)。
    hubble_km_s_kpc = 0.1 * cosmology["h"] * np.sqrt(
        cosmology["omega_m0"] * (1.0 + redshift_grid) ** 3
        + cosmology["omega_lambda0"]
    )

    # 临界密度 rho_crit(z)=3 H(z)^2/(8 pi G)。当前 H 和 G 的单位组合
    # 会让 rho_crit 直接得到 M_sun/kpc^3，不需要额外长度或质量换算。
    critical_density_msun_kpc3 = (
        3.0
        * hubble_km_s_kpc**2
        / (8.0 * np.pi * GRAVITATIONAL_CONSTANT_KPC_KM2_S2_MSUN)
    )

    # M_200c=(4 pi/3)*200*rho_crit*R_vir^3，反解得到 R_vir；
    # 浓度定义 C=R_vir/R_s，因此 R_s=R_vir/C。
    virial_radius_kpc = np.cbrt(
        3.0
        * mass_tracks_msun
        / (4.0 * np.pi * VIRIAL_OVERDENSITY * critical_density_msun_kpc3)
    )
    scale_radius_kpc = virial_radius_kpc / concentration_tracks

    return (
        mass_tracks_msun,
        concentration_tracks,
        virial_radius_kpc,
        scale_radius_kpc,
    )


def main():
    """生成 Fig. 9 及其曲线数据。"""
    make_figure9(Path(__file__).resolve().parent)
    print("Created fig9_reproduction.png")
    print("Created fig9_reproduction.csv")


if __name__ == "__main__":
    main()
