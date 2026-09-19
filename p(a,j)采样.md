# Appendix D 的 $P(a,j)$ 与联合采样程序

## 1. 这份文档解决什么问题

目标论文 Appendix D 给出了原初黑洞双星形成时的半长轴 $a$、重标度角动量
$j$ 和偏心率 $e$ 的分布。程序实现位于：

    p_a_j.py

这份文档分两层说明：

1. Appendix D 的 D1-D5 各公式分别描述什么；
2. 怎样把连续联合密度 $P(a,j)$ 变成有限数量的相关随机样本。

第二点最重要。$a$ 和 $j$ 不是两个互相独立的随机变量，因为 Eq. (D3)
中的角动量尺度依赖 $x/\bar x$，而 $x$ 又通过 Eq. (D5) 依赖 $a$。
因此不能先从某个固定的 $P(a)$ 采样 $a$，再从一个与 $a$ 无关的
$P(j)$ 采样 $j$。

当前程序采用的联合采样顺序为

$$
P(a,j)
\longrightarrow
p_a(a)=\int_0^1P(a,j)\,dj
\longrightarrow
a=F_a^{-1}(u_a)
\longrightarrow
j=F_{j|a}^{-1}(u_j)
\longrightarrow
e=\sqrt{1-j^2}.
$$

这里 $u_a,u_j$ 是两个独立的 $[0,1)$ 均匀随机数，但最后得到的 $a,j$
并不独立，因为 $j$ 的条件分布依赖已经抽到的 $a$。

---

## 2. 变量、单位和默认参数

### 2.1 轨道变量

- $a$：PBH 双星形成时的半长轴，程序统一使用 pc；
- $e$：轨道偏心率，满足 $0\leq e<1$；
- $j$：重标度角动量，

$$
j=\sqrt{1-e^2},
$$

满足 $0<j<1$。

这里的 $j$ 不是带有物理单位的轨道角动量，而是一个无量纲变量。

当 $j\rightarrow0$ 时，

$$
e\rightarrow1,
$$

对应高度偏心的轨道。当 $j\rightarrow1$ 时，

$$
e\rightarrow0,
$$

对应近圆轨道。

### 2.2 PBH 和宇宙学变量

- $m_{\rm PBH}$：单个 PBH 的质量，单位 $M_\odot$；
- $f_{\rm PBH}$：暗物质中 PBH 所占比例；
- $\rho_{\rm eq}$：物质-辐射相等时期的平均物质质量密度；
- $\sigma_{\rm eq}\simeq0.005$：相等时期大尺度高斯密度涨落的方差参数；
- $\alpha=0.1$：半长轴与初始 PBH 间距关系中的系数；
- $z_{\rm eq}\simeq3450$：物质-辐射相等红移。

脚本的默认参数为：

    DEFAULT_OMEGA_M0 = 0.308
    DEFAULT_H = 0.678
    DEFAULT_Z_EQ = 3450.0
    DEFAULT_SIGMA_EQ = 0.005
    DEFAULT_ALPHA = 0.1

其中 $\Omega_{m0}$ 和 $h$ 与项目中 Ludlow16 使用的 Planck 参数一致。

需要特别注意：论文 Appendix D 只给出 $z_{\rm eq}\simeq3450$，没有单独说明
计算 $\rho_{\rm eq}$ 时使用的 $\Omega_{m0}$ 和 $h$。因此
p_a_j.py 允许调用者直接传入自己的 $\rho_{\rm eq}$，默认值只是当前项目明确
采用的一种口径。

---

## 3. 从今天的宇宙学参数得到 $\rho_{\rm eq}$

函数 critical_density_z0_msun_pc3 先计算今天的临界密度：

$$
\rho_{\rm crit,0}
=
\frac{3H_0^2}{8\pi G},
\qquad
H_0=100h\ {\rm km\,s^{-1}\,Mpc^{-1}}.
$$

脚本使用

$$
G=4.30091727\times10^{-3}
\ {\rm pc}\,({\rm km\,s^{-1}})^2M_\odot^{-1},
$$

所以输出可以直接写成 $M_\odot\,{\rm pc}^{-3}$。

函数 matter_density_at_equality_msun_pc3 再使用非相对论物质的红移关系：

$$
\rho_m(z)
=
\Omega_{m0}\rho_{\rm crit,0}(1+z)^3.
$$

因此

$$
\rho_{\rm eq}
=
\Omega_{m0}\rho_{\rm crit,0}(1+z_{\rm eq})^3.
$$

默认参数下程序得到

$$
\rho_{\rm eq}\simeq1614.97\ M_\odot\,{\rm pc}^{-3}.
$$

论文称 $\rho_{\rm eq}$ 为能量密度；当前代码使用质量密度单位，这是在
$c=1$ 的公式与天体物理质量单位之间进行的明确换算。

---

## 4. Eq. (D1)：平均 PBH 间距

论文定义物质-辐射相等时期的平均 PBH 间距为

$$
\bar x
=
\left(
\frac{3m_{\rm PBH}}
{4\pi f_{\rm PBH}\rho_{\rm eq}}
\right)^{1/3}.
\tag{D1}
$$

对应函数为：

    mean_pbh_separation_pc(...)

默认的

$$
m_{\rm PBH}=30M_\odot,\qquad f_{\rm PBH}=1
$$

给出

$$
\bar x\simeq0.1643\ {\rm pc}.
$$

从公式可以直接看出：

$$
\bar x\propto m_{\rm PBH}^{1/3},
\qquad
\bar x\propto f_{\rm PBH}^{-1/3},
\qquad
\bar x\propto\rho_{\rm eq}^{-1/3}.
$$

PBH 越重，平均间距越大；PBH 占比或背景密度越高，平均间距越小。

---

## 5. Eq. (D5)：由半长轴反推初始间距

论文给出

$$
x(a)
=
\left(
\frac{3am_{\rm PBH}}
{4\pi\alpha\rho_{\rm eq}}
\right)^{1/4}.
\tag{D5}
$$

对应函数为：

    separation_from_semimajor_axis_pc(...)

它把双星形成前的 PBH 空间间距 $x$ 与形成后的半长轴 $a$ 联系起来。

将 D1 和 D5 联立，可以得到一个非常有用的无量纲变量：

$$
X(a)
\equiv
\left(\frac{x(a)}{\bar x}\right)^3
=
\left(
\frac{af_{\rm PBH}}{\alpha\bar x}
\right)^{3/4}.
$$

于是

$$
a
=
\frac{\alpha\bar x}{f_{\rm PBH}}X^{4/3}.
$$

这个关系解释了 D4 中指数项

$$
\exp[-X(a)].
$$

它来自初始 PBH 泊松空间分布中最近邻距离的指数概率。

---

## 6. Eqs. (D2)-(D3)：给定 $a$ 时的角动量分布

论文写出

$$
P(j)
=
\frac{y(j)^2}
{j[1+y(j)^2]^{3/2}},
\tag{D2}
$$

其中

$$
y(j)
=
\frac{j}
{
0.5
\sqrt{1+\sigma_{\rm eq}^2/f_{\rm PBH}^2}
(x/\bar x)^3
}.
\tag{D3}
$$

为了让公式和代码更清楚，定义角动量尺度

$$
j_0(a)
=
0.5
\sqrt{1+\frac{\sigma_{\rm eq}^2}{f_{\rm PBH}^2}}
\left(\frac{x(a)}{\bar x}\right)^3.
$$

于是

$$
y(j)=\frac{j}{j_0(a)},
$$

并且

$$
P(j\mid a)
=
\frac{(j/j_0)^2}
{j[1+(j/j_0)^2]^{3/2}}.
$$

程序中的函数对应关系为：

- angular_momentum_scale：计算 $j_0(a)$；
- angular_momentum_pdf：计算 D2-D3；
- angular_momentum_cdf：计算 D2 的解析积分。

### 6.1 为什么 $P(j)$ 实际依赖 $a$

D2 表面上写成 $P(j)$，但 D3 中含有 $x/\bar x$。根据 D5，

$$
x=x(a),
$$

因此更完整的写法应当是

$$
P(j\mid a).
$$

这正是 $a$ 与 $j$ 相关的来源。不同半长轴对应不同的 $j_0$，所以不能
为所有 $a$ 使用同一条固定的 $P(j)$。

### 6.2 D2 的解析 CDF

将

$$
y=\frac{j}{j_0}
$$

代入 D2，可以得到从 0 到 $j$ 的积分：

$$
F_{\rm raw}(j\mid a)
=
\int_0^jP(j'\mid a)\,dj'
=
1-\frac{1}{\sqrt{1+(j/j_0)^2}}.
$$

这就是 angular_momentum_cdf 的公式。

若积分上限取无穷大，则

$$
\lim_{j\rightarrow\infty}F_{\rm raw}(j\mid a)=1.
$$

但是论文的物理变量要求

$$
0<j<1.
$$

因此 D2 的原始密度在有限区间内的积分一般不是 1，而是

$$
Q(a)
\equiv
F_{\rm raw}(1\mid a)
=
1-\frac{1}{\sqrt{1+j_0(a)^{-2}}}.
$$

真正用于 $0<j<1$ 条件采样的归一化密度为

$$
p(j\mid a,0<j<1)
=
\frac{P(j\mid a)}{Q(a)}.
$$

angular_momentum_pdf 的 normalize_on_unit_interval 参数就是用来区分：

- False：返回论文印刷的 D2 原始密度；
- True：返回在 $0<j<1$ 上重新归一化的条件密度。

---

## 7. Eq. (D4)：半长轴和角动量的联合密度

论文给出

$$
P(a,j)
=
\frac{3a^{-1/4}}{4}
\left(
\frac{f_{\rm PBH}}{\alpha\bar x}
\right)^{3/4}
P(j\mid a)
\exp\left[
-\left(\frac{x(a)}{\bar x}\right)^3
\right].
\tag{D4}
$$

程序函数为：

    joint_orbital_pdf(...)

为了突出结构，可以把 D4 写成

$$
P(a,j)=B(a)P(j\mid a),
$$

其中

$$
B(a)
=
\frac{3a^{-1/4}}{4}
\left(
\frac{f_{\rm PBH}}{\alpha\bar x}
\right)^{3/4}
\exp[-X(a)].
$$

$B(a)$ 的单位为 ${\rm pc}^{-1}$，而 $P(j\mid a)$ 对无量纲 $j$
而言是无量纲密度，所以联合密度对 $a$ 的单位也是 ${\rm pc}^{-1}$。

### 7.1 joint_orbital_pdf 为什么不自动归一化

joint_orbital_pdf 忠实返回论文印刷的 D4。它没有假装自己已经在

$$
10^{-6}<a<1\ {\rm pc},
\qquad
0<j<1
$$

这个有限矩形区域内归一化。

原因是“原始理论密度”和“为了绘图而截断后的采样密度”是两件不同的事情。
如果在基础公式函数里悄悄做有限区间归一化，以后改变 $a$ 范围时就很容易
不知道归一化发生在哪里。

---

## 8. 从联合密度得到 $a$ 的边缘密度

联合采样的第一步不是随便猜一条 $P(a)$，而是把 D4 对允许的 $j$ 积分：

$$
p_a^{\rm raw}(a)
=
\int_0^1P(a,j)\,dj.
$$

由于

$$
P(a,j)=B(a)P(j\mid a)
$$

且 D2 的积分已经有解析式，所以

$$
p_a^{\rm raw}(a)
=
B(a)Q(a),
$$

其中

$$
Q(a)=F_{\rm raw}(1\mid a).
$$

这就是 semi_major_axis_marginal_pdf 的实现：

    return semi_major_axis_factor * angular_momentum_cdf(
        1.0,
        j_scale,
    )

这个 $Q(a)$ 因子很重要。若把 D2 当成已经在 $0<j<1$ 上归一化而直接丢掉
$Q(a)$，得到的是另一种半长轴分布。

在论文使用的有限半长轴范围

$$
a_{\min}=10^{-6}\ {\rm pc},
\qquad
a_{\max}=1\ {\rm pc}
$$

内，还需要再计算整体归一化：

$$
N_a
=
\int_{a_{\min}}^{a_{\max}}
p_a^{\rm raw}(a)\,da.
$$

最终用于抽取 $a$ 的概率密度是

$$
p_a(a)
=
\frac{p_a^{\rm raw}(a)}{N_a}.
$$

---

## 9. 半长轴的数值逆 CDF 采样

sample_initial_orbital_parameters 的第一部分负责采样 $a$。

### 9.1 为什么使用对数网格

论文的范围跨越六个数量级：

$$
10^{-6}\ {\rm pc}<a<1\ {\rm pc}.
$$

若使用线性网格，小 $a$ 区域只有很少的网格点，无法分辨概率密度快速变化的
部分。因此程序建立

    semi_major_axis_grid_pc = np.geomspace(
        minimum_semi_major_axis_pc,
        maximum_semi_major_axis_pc,
        semi_major_axis_grid_point_count,
    )

默认使用 16384 个对数等距节点。

### 9.2 对数网格不等于对 $d\ln a$ 积分

虽然节点是对数分布的，但目标概率仍然是

$$
p_a(a)\,da.
$$

因此代码调用：

    semi_major_axis_cdf = cumulative_trapezoid(
        marginal_density,
        x=semi_major_axis_grid_pc,
        initial=0.0,
    )

这里显式传入

    x=semi_major_axis_grid_pc

表示梯形积分使用真实的 $da$ 间隔。

如果改为对 $\ln a$ 积分，就必须使用

$$
da=a\,d\ln a
$$

并在被积函数中额外乘上 Jacobian $a$。当前程序没有遗漏这个 Jacobian；
它选择的是更直接的“在非均匀 $a$ 节点上对 $da$ 积分”。

### 9.3 构造并归一化数值 CDF

梯形累计积分得到

$$
C_i
\simeq
\int_{a_{\min}}^{a_i}
p_a^{\rm raw}(a)\,da.
$$

随后执行

    semi_major_axis_cdf /= semi_major_axis_cdf[-1]

使最后一个网格点满足

$$
F_a(a_{\max})=1.
$$

### 9.4 用 np.interp 实现逆 CDF

程序生成 $N$ 个均匀随机数：

    u_a = rng.random(sample_count)

然后使用

    semi_major_axis_pc = np.interp(
        u_a,
        semi_major_axis_cdf,
        semi_major_axis_grid_pc,
    )

完成

$$
a=F_a^{-1}(u_a).
$$

np.interp 的第一个已知坐标是 CDF，第二个已知坐标是半长轴网格，所以这里做的
不是普通的“由 $a$ 求概率”，而是反过来“由累计概率求 $a$”。

---

## 10. 给定每个 $a$ 后采样 $j$

得到每一个 $a_i$ 后，程序重新计算：

$$
x_i=x(a_i),
\qquad
j_{0,i}
=
0.5
\sqrt{1+\sigma_{\rm eq}^2/f_{\rm PBH}^2}
\left(\frac{x_i}{\bar x}\right)^3.
$$

因此不同样本有不同的角动量尺度 $j_{0,i}$。

### 10.1 截断后的目标 CDF

对于每个 $a_i$，先计算

$$
Q_i
=
F_{\rm raw}(1\mid a_i).
$$

再生成均匀随机数 $u_{j,i}\in[0,1)$，令

$$
T_i=u_{j,i}Q_i.
$$

代码对应为：

    probability_below_one = angular_momentum_cdf(1.0, j_scale)
    target_cdf = rng.random(sample_count) * probability_below_one

这样 $T_i$ 均匀覆盖原始 CDF 的

$$
[0,F_{\rm raw}(1\mid a_i))
$$

区间，等价于从已经归一化的 $0<j<1$ 条件分布中抽样。

### 10.2 D2 的解析逆 CDF

从

$$
T
=
1-\frac{1}{\sqrt{1+y^2}}
$$

开始，移项得到

$$
\frac{1}{\sqrt{1+y^2}}=1-T.
$$

两边平方并求解：

$$
1+y^2=(1-T)^{-2},
$$

$$
y
=
\sqrt{(1-T)^{-2}-1}.
$$

因为

$$
j=j_0y,
$$

所以最终的逆变换为

$$
j
=
j_0
\sqrt{(1-T)^{-2}-1}.
$$

程序代码正是：

    inverse_cdf_denominator = 1.0 - target_cdf
    y = np.sqrt(inverse_cdf_denominator ** (-2.0) - 1.0)
    angular_momentum = j_scale * y

最后使用

    np.nextafter(1.0, 0.0)

限制极少数浮点舍入产生的 $j=1$。这不是物理截断参数，而是防止下一步计算
$\sqrt{1-j^2}$ 时因为浮点误差越过定义域。

---

## 11. 从 $j$ 转换到偏心率 $e$

论文定义

$$
j=\sqrt{1-e^2}.
$$

因此

$$
e=\sqrt{1-j^2}.
$$

程序提供两个方向的函数：

    angular_momentum_from_eccentricity(eccentricity)
    eccentricity_from_angular_momentum(j)

采样器最终调用：

    eccentricity = eccentricity_from_angular_momentum(
        angular_momentum
    )

需要强调：不能把 $j$ 的概率密度直接当成 $e$ 的概率密度。如果要解析地
由 $P(j)$ 变换成 $P(e)$，必须包含 Jacobian：

$$
P(e\mid a)
=
P(j(e)\mid a)
\left|\frac{dj}{de}\right|,
$$

其中

$$
\left|\frac{dj}{de}\right|
=
\frac{e}{\sqrt{1-e^2}}
=
\frac{e}{j}.
$$

当前程序先正确采样 $j$，再逐样本确定性地转换为 $e$，因此 Jacobian 已由
随机变量变换自动体现，不需要在样本上再乘一次。

---

## 12. 为什么这种方法保留了 $a$-$j$ 相关性

完整算法可以写成：

1. 由 D1 求 $\bar x$；
2. 在 $a$ 网格上由 D5 求 $x(a)$；
3. 由 D3 求每个 $a$ 对应的 $j_0(a)$；
4. 把 D4 对 $0<j<1$ 积分，得到 $p_a(a)=B(a)Q(a)$；
5. 从 $p_a(a)$ 的数值 CDF 抽取 $a_i$；
6. 对每个 $a_i$ 单独计算 $j_{0,i}$；
7. 从 $p(j\mid a_i,0<j<1)$ 抽取 $j_i$；
8. 计算 $e_i=\sqrt{1-j_i^2}$。

虽然步骤 5 和步骤 7 使用两个独立均匀随机数，但 $j_i$ 的变换函数含有
$a_i$，所以联合样本仍然服从相关分布。

错误的独立采样方式会写成：

    a = sample_from_one_fixed_p_a(...)
    j = sample_from_one_fixed_p_j(...)

它相当于假设

$$
P(a,j)=P(a)P(j),
$$

而 D3-D5 明确说明这个分解一般不成立。

---

## 13. 程序的返回对象

sample_initial_orbital_parameters 返回 OrbitalSamples：

    @dataclass(frozen=True)
    class OrbitalSamples:
        semi_major_axis_pc: np.ndarray
        angular_momentum: np.ndarray
        eccentricity: np.ndarray

三个数组的相同下标表示同一颗 PBH 双星。例如：

    samples.semi_major_axis_pc[i]
    samples.angular_momentum[i]
    samples.eccentricity[i]

必须作为同一个三元组理解，不能分别排序后再组合，否则会破坏 $a$-$j$
相关性。

frozen=True 表示数据类的字段引用创建后不能被重新赋值。数组内容在 NumPy
层面仍可修改，因此分析阶段仍应避免原地打乱单个字段。

---

## 14. 面向后续模拟的类接口

AppendixDOrbitalDistribution 保存一整套分布参数：

    distribution = AppendixDOrbitalDistribution(
        pbh_mass_msun=30.0,
        f_pbh=1.0,
        sigma_eq=0.005,
        alpha=0.1,
        minimum_semi_major_axis_pc=1.0e-6,
        maximum_semi_major_axis_pc=1.0,
    )

然后通过统一接口采样：

    rng = np.random.default_rng(12345)
    samples = distribution.sample(100_000, rng)

这符合 PLAN.md 约定的

    OrbitalDistribution.sample(n, rng)

形式。将来如果换成其他 $P(a,e)$，轨道演化模块只需要接收同样结构的样本，
不需要知道 Appendix D 的具体逆 CDF 实现。

---

## 15. 为什么随机数生成器必须从外部传入

程序没有在函数内部调用

    np.random.seed(...)

也没有使用旧式全局随机函数。调用者必须明确创建：

    rng = np.random.default_rng(20260914)

这样做有三个优点：

1. 固定种子时可以精确重复同一组样本；
2. 不会悄悄改变项目其他模块的全局随机状态；
3. 多个晕、壳层或并行任务可以拥有彼此独立的随机数流。

例如：

    distribution = AppendixDOrbitalDistribution()

    first = distribution.sample(
        10_000,
        np.random.default_rng(14),
    )
    second = distribution.sample(
        10_000,
        np.random.default_rng(14),
    )

first 和 second 会得到相同样本。如果连续两次调用同一个 rng 对象：

    rng = np.random.default_rng(14)
    first = distribution.sample(10_000, rng)
    second = distribution.sample(10_000, rng)

第二次调用会从同一个随机数流的后续状态继续，因此样本不会重复。这通常是模拟
多个批次时希望得到的行为。

---

## 16. 一次完整调用示例

    import numpy as np

    import p_a_j as paj

    distribution = paj.AppendixDOrbitalDistribution(
        pbh_mass_msun=30.0,
        f_pbh=1.0,
        minimum_semi_major_axis_pc=1.0e-6,
        maximum_semi_major_axis_pc=1.0,
    )

    rng = np.random.default_rng(20260914)
    samples = distribution.sample(100_000, rng)

    a_pc = samples.semi_major_axis_pc
    j = samples.angular_momentum
    e = samples.eccentricity

    a_au = a_pc * paj.AU_PER_PC

这里的 a_pc、j、e 和 a_au 都是一维数组；第 $i$ 个位置始终属于同一颗双星。

若只想评价理论联合密度而不采样，可以使用：

    density = paj.joint_orbital_pdf(
        semi_major_axis_pc=a_values,
        j=j_values,
        pbh_mass_msun=30.0,
        f_pbh=1.0,
    )

---

## 17. 论文 Fig. 19 中仍存在的实现歧义

### 17.1 样本数矛盾

正文说逆变换采样使用

$$
N_{\rm sample}=10^4,
$$

但 Fig. 19 图注写

$$
N=10^5.
$$

当前脚本不擅自选择其中一个；sample_count 由绘图 Notebook 明确传入。若目标是
复现图中的计数纵轴，优先尝试 $10^5$，同时在图注中记录正文的 $10^4$
说法。

### 17.2 “先采样 $j$”时缺少 $x/\bar x$

正文说先从 $P(j)$ 采样 $j$，再由联合分布采样 $a$。但是 D3 中的
$P(j)$ 需要 $x/\bar x$，而 D5 又让 $x$ 依赖 $a$。在尚未抽到 $a$
时，作者没有说明这里的 $x/\bar x$ 取什么值。

当前程序没有猜测一个固定 $x/\bar x$，而是从完整联合密度出发，先采样正确的
$a$ 边缘，再采样 $j\mid a$。

### 17.3 有限区间归一化顺序

论文给出

$$
10^{-6}<a<1\ {\rm pc},
\qquad
0<j<1,
$$

但没有说明是：

1. 先使用无限区间理论分布采样，再拒绝超出范围的样本；
2. 还是先把联合分布限制到有限矩形，再整体归一化。

在拒绝采样严格实现且接受条件相同时，两种方法应给出相同的条件分布；但如果作者
分两步分别归一化 $P(j)$ 和 $P(a,j)$，结果可能不同。

### 17.4 平均半长轴差异

论文声称

$$
\langle a\rangle\simeq7287\ {\rm AU}.
$$

当前代码采用：

- Planck/Ludlow16 的 $\Omega_{m0}=0.308,\ h=0.678$；
- $z_{\rm eq}=3450$；
- $m_{\rm PBH}=30M_\odot$；
- $f_{\rm PBH}=1$；
- D4 在给定有限矩形内的严格联合截断归一化。

使用 $10^5$ 个固定种子样本时得到的平均值约为

$$
\langle a\rangle\simeq2.14\times10^3\ {\rm AU}.
$$

这说明论文 Fig. 19 的实际代码还可能使用了不同的 $\rho_{\rm eq}$ 口径、不同的
截断顺序，或正文没有说明的采样步骤。当前不能把 7287 AU 当作已经由 D1-D5
和已公开参数唯一推出的结果。

因此代码保留物理参数接口和实现说明，没有通过修改 $\rho_{\rm eq}$、$\alpha$
或额外乘常数来强行得到 7287 AU。

---

## 18. 常见错误

### 错误 1：独立采样 $a$ 和 $j$

D3 中的 $j_0$ 依赖 $a$，所以独立采样会破坏联合分布。

### 错误 2：认为 D2 自动在 $0<j<1$ 上归一化

D2 在 $0<j<\infty$ 上归一化。限制到 $j<1$ 后必须除以 $Q(a)$；同时
$a$ 的边缘分布中必须保留 $Q(a)$。

### 错误 3：使用对数网格后直接对数组下标求和

对数网格相邻点的 $\Delta a$ 不相等。必须使用真实 a_grid 作为积分横坐标，
或者改写成对 $d\ln a$ 积分并显式乘 Jacobian $a$。

### 错误 4：把 $j$ 均匀分布当作偏心率均匀分布

$$
e=\sqrt{1-j^2}
$$

是非线性变换。即使 $j$ 均匀，$e$ 也不会均匀；而 Appendix D 中的 $j$
本身也不均匀。

### 错误 5：转换 $P(j)$ 到 $P(e)$ 时重复乘 Jacobian

若直接变换解析密度，需要 Jacobian。若先抽取 $j$ 样本再逐个计算 $e$，变量
变换已经自动发生，不能再给样本附加一次 Jacobian 权重。

### 错误 6：分别排序 $a$、$j$、$e$

三个数组相同下标代表同一颗双星。分别排序会破坏轨道参数之间的物理对应关系。

### 错误 7：把 $\rho_{\rm eq}$ 的不确定口径隐藏在常数中

论文没有公开完整数值实现。若为了匹配 Fig. 19 调整 $\rho_{\rm eq}$，必须把
所用定义和数值写入绘图 Notebook，不能悄悄修改模块默认值。

---

## 19. 与 Fig. 19 绘图 Notebook 的职责划分

p_a_j.py 只负责：

- 公式 D1-D5；
- 有限区间联合采样；
- 随机数接口；
- $a,j,e$ 的返回结构。

后续 Fig. 19 Notebook 应负责：

- 明确选择 $N=10^4$ 还是 $10^5$；
- 选择直方图 bin 数与范围；
- 对半长轴和偏心率分别画计数直方图；
- 使用对数纵轴；
- 输出平均半长轴和所用 $\rho_{\rm eq}$；
- 将程序结果与论文原始矢量图比较；
- 把 7287 AU 的差异保留为待排查项。

这种分工可以保证绘图样式或样本数变化不会修改基础概率分布模块。

---

## 20. 最简调用链总结

从程序入口看，推荐使用：

    distribution = AppendixDOrbitalDistribution(...)
    samples = distribution.sample(sample_count, rng)

内部真实调用关系为：

$$
\rho_{\rm eq}
\rightarrow
\bar x
\rightarrow
x(a)
\rightarrow
j_0(a)
\rightarrow
p_a(a)
\rightarrow
F_a(a)
\rightarrow
a_i
\rightarrow
F_{j|a_i}(j)
\rightarrow
j_i
\rightarrow
e_i.
$$

理解这条调用链后，Fig. 19 的两个直方图就不再是两个互不相关的随机分布：
上图显示联合样本在 $a$ 方向的投影，下图显示同一批联合样本经
$e=\sqrt{1-j^2}$ 变换后的投影。

