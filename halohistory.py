"""暗物质晕质量吸积史：论文 Fig. 9 的逐步复现脚本。

这个文件对应 9 月 9 日的学习任务。它不重复实现 Ludlow16 浓度公式，
而是导入 ``haloconcentration.py`` 中已经完成并核对过的函数。

当前阶段完成两件事：

1. 为 Fig. 9 的一组现今晕质量计算 M(z)；
2. 使用该数据矩阵绘制并保存 Fig. 9。

后续阶段才会加入 M(z) -> C(z) -> R_vir(z) -> R_s(z) 的数据流。
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


def calculate_mass_accretion_tracks(
    z,
    present_day_masses_msun=FIG9_PRESENT_DAY_MASSES_MSUN,
    preserve_mass_rank=False,
):
    """计算 Fig. 9 中不同现今质量晕的平均质量吸积轨迹。

    Parameters
    ----------
    z : float 或 array-like
        需要回溯的红移。z=0 表示今天；z 越大，表示回溯到越早的宇宙。
    present_day_masses_msun : array-like
        每条轨迹在 z=0 时的晕质量 M0，单位统一为太阳质量 M_sun。
    preserve_mass_rank : bool
        若为 True，在每个红移按 M0 的质量秩重新排列 Eq. (31) 的结果，使较大
        M0 对应的代表性晕质量不小于较小 M0。该操作是 Fig. 9 的族群秩序约束，
        不是 Ludlow16 或 Correa15 原始公式；False 返回未经处理的公式轨迹。

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

    mass_tracks_msun = hc.mass_accretion_history(
        present_day_masses_msun[:, None], z[None, :], model="ludlow16"
    )
    if not preserve_mass_rank:
        return mass_tracks_msun

    # 独立的中位 MAH 公式在高质量、高红移外推区会发生秩交换。Fig. 9 若解释为
    # 按现今质量排序的代表性晕族群，应在每个红移保留该质量秩；原始值仍可通过
    # preserve_mass_rank=False 完整取得。
    mass_order = np.argsort(present_day_masses_msun)
    ordered_tracks = np.empty_like(mass_tracks_msun)
    ordered_tracks[mass_order] = np.sort(mass_tracks_msun, axis=0)
    return ordered_tracks


def make_figure9(output_directory):
    """使用质量吸积轨迹矩阵绘制并保存论文 Fig. 9。"""
    output_directory = Path(output_directory)
    z = np.geomspace(0.1, 12.0, 400)
    mass_tracks_msun = calculate_mass_accretion_tracks(
        z, preserve_mass_rank=True
    )

    figure, axis = plt.subplots(figsize=(7.0, 4.8))
    line_styles = ("-", "--", "-.", ":")
    for index, mass0_msun in enumerate(FIG9_PRESENT_DAY_MASSES_MSUN):
        exponent = int(np.log10(mass0_msun))
        axis.plot(
            z,
            mass_tracks_msun[index],
            linestyle=line_styles[index % len(line_styles)],
            label=rf"$M=10^{{{exponent}}}M_\odot$",
        )

    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(0.1, 12.0)
    axis.set_xlabel(r"Redshift $z$")
    axis.set_ylabel(r"$M(z)\;[M_\odot]$")
    axis.legend(loc="upper right", fontsize="small")
    figure.tight_layout()
    figure.savefig(output_directory / "fig9_reproduction(1).png", dpi=220)
    plt.close(figure)
    return z, mass_tracks_msun


# TODO（再下一步）：加入 R_vir(z) 与 R_s(z)，建立完整晕结构数据流。


def main():
    """生成 Fig. 9。"""
    make_figure9(Path(__file__).resolve().parent)
    print("Created fig9_reproduction(1).png")


if __name__ == "__main__":
    main()
