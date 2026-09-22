# 第 IV 节数值模拟：从输入、轨道演化到并合率输出

## 0. 本文档的范围

本文档整理目标论文第 IV 节及其直接相关的附录 D、E，说明其中“原初黑洞双星—单体相互作用”数值模拟的完整逻辑。重点不是逐句翻译，而是把论文散落在正文、公式、图注和表格中的设置重组成一条可实现的数据流。

为避免把复现代码中的选择误写成论文事实，下文使用三类标记：

- **论文明确给出**：可以直接由正文、公式或表格确认；
- **由论文设定直接推出**：是实现时必需、且可由论文定义推导的量，但论文未单独写成公式；
- **论文未说明/存在歧义**：不能在实现时静默补齐，必须作为项目自己的数值选择记录。

特别说明：论文第 IV 节采用 **Ludlow16 浓度关系**。本项目后续默认采用 **Prada12-HMF 吸积史/浓度接口**，这是复现工程的可替换模型选择，不应倒写成论文原始设置。

---

## 1. 一句话说明模拟的核心

这不是逐个解析 PBH 散射过程的直接 $N$-体模拟，而是一个**分层的 Monte Carlo 转移估计器**：

1. 用质量吸积史和 NFW 模型给出每个红移、每个壳层的平滑环境；
2. 从早期 PBH 双星分布中抽取相同数量的轨道样本；
3. 只保留在当地环境中属于硬双星的样本；
4. 用平均化的双星—单体硬化公式与引力波反作用公式演化 $a,e$；
5. 统计该环境中样本在一个时间片内并合的比例；
6. 再乘以该壳层真实拥有的双星数，得到物理并合率；
7. 最后用 halo mass function 对不同 halo 质量积分，得到共动体积并合率。

其最核心的数学结构是

$$
\boxed{
\text{某环境中的并合概率}
\times
\text{该环境中的真实双星数}
=
\text{该环境贡献的并合数}
}
$$

因此，“每个壳层都模拟同样多的双星”只是为了控制 Monte Carlo 误差，绝不表示各壳层真实双星数相同。

---

## 2. 模拟对象与没有模拟的对象

### 2.1 被模拟的对象

论文跟踪单个 PBH 双星的两个轨道变量：

$$
X=(a,e),
$$

其中 $a$ 是半长轴，$e$ 是偏心率。每个样本还附带以下环境标签：

$$
(M_0,t_k,z_k,i),
$$

分别表示该 halo 的今天质量、当前时间/红移以及所在壳层编号。

论文生产计算采用单色 PBH 质量：

$$
m_1=m_2=m=30\,M_\odot.
$$

等质量假设使交换相互作用不需要单独处理。论文也把 binary-binary 相互作用近似为“较硬的双星与最近一个 PBH”的 binary-single 相互作用。

### 2.2 没有被直接模拟的内容

以下过程没有通过逐次散射或 $N$-体动力学显式求解：

- 每一次三体接近的冲量、入射参数和散射角；
- halo 内 PBH 的离散相空间分布；
- 壳层之间的径向迁移；
- 双星被踢出 halo 后的继续演化；
- 双星总体因并合、解体或逃逸造成的自洽耗尽；
- binary-binary 的四体动力学。

这些效应被平滑的 NFW 背景、经验硬化系数 $H,K$ 和若干近似取代。论文明确忽略 PBH 逃逸，因此低质量 halo 的并合率可能被高估。

---

## 3. 输入层：开始模拟前必须具备什么

模拟输入可分为五层。

| 输入层 | 主要变量 | 作用 |
|---|---|---|
| 宇宙学 | $H_0,\Omega_m,\rho_{\rm crit}(z)$，时间—红移关系 | 把红移轨迹转换成可积分的时间网格，并定义 virial 量 |
| halo 演化 | $M_0,M(z),C(M,z),R_{\rm vir}(z)$ | 建立随时间变化的 NFW halo |
| 空间离散 | $N_{\rm shell},R_i(t)$ | 将 halo 分解成若干固定编号、边界随时间变化的层 |
| PBH 总体 | $m,f_{\rm PBH},f_{\rm binary},f_{\rm single}$ | 决定环境密度和各壳层真实双星数 |
| 双星初态与积分参数 | $P(a,j),N_{\rm sample},\Delta t,\delta t$ | 生成 Monte Carlo 样本并演化轨道 |

### 3.1 halo 质量吸积史

论文使用

$$
M(z)=M_0(1+z)^\alpha e^{\beta z},
\tag{31}
$$

其中 $M_0=M(z=0)$，而 $\alpha,\beta$ 由附录 C 的吸积史拟合确定。逻辑顺序必须是

$$
M_0
\longrightarrow M(z)
\longrightarrow C[M(z),z]
\longrightarrow R_{\rm vir}(z),r_s(z),\rho_s(z).
$$

不能把 $z=0$ 的 halo 结构参数直接用于所有红移。

论文在这里采用 Ludlow16 浓度模型。工程上可以把浓度/吸积史做成通用接口，但替换为 Prada12-HMF 后得到的是“项目基线结果”，不是严格的论文原设结果。

### 3.2 起始红移

论文脚注给出的规则是：

- 今天质量 $M_0\gtrsim5\times10^4M_\odot$ 的 halo 从 $z=12$ 开始；
- 更小的 halo 从较低红移开始，以保证任一时刻 halo 中至少约有 30 个 PBH。

从 $z=12$ 到 $z=0$ 约对应 13.4 Gyr。这里“至少约 30 个 PBH”是低质量端的数值适用性限制，不是新的物理形成阈值。

### 3.3 PBH 比例

第 IV 节的基准计算取 $f_{\rm PBH}=1$，并设置

$$
f_{\rm binary}=f_{\rm single}=0.5,
\qquad
f_{\rm binary}+f_{\rm single}=1.
$$

论文把双星 PBH 与单体 PBH 的质量密度相加作为散射环境：

$$
\rho_{\rm env}
=\rho_{\rm PBH,binary}+\rho_{\rm PBH,single}
=\rho_{\rm NFW}.
\tag{21}
$$

宽双星在环境密度中被当成两个独立 PBH。这里的 $f_{\rm binary}$ 是论文的总体比例设定；从质量比例换算成“真实双星系统数”时还必须除以每个系统的总质量 $2m$。

---

## 4. halo 的层级划分与局部环境

### 4.1 壳层数量

论文采用分质量档的空间分辨率：

| halo 质量档 | 空间划分 |
|---|---:|
| $mathcal O(10^3M_\odot)$ | 1 个球体 |
| $mathcal O(10^4M_\odot)$ | 1 个内球 + 1 个外壳，共 2 层 |
| $mathcal O(10^5M_\odot)$ | 1 个内球 + 2 个外壳，共 3 层 |
| $mathcal O(10^6M_\odot)$ | 5 层 |
| $M\gtrsim10^7M_\odot$ | 10 层 |

论文正文没有一句话直接写“壳层数由 $z=0$ 的质量 $M_0$ 决定”。不过附录 E 的表 III 中，$M_0=10^6M_\odot$ 的轨迹在 $z=8$ 已降到 $M(z)=7.3\times10^4M_\odot$，表中仍然沿用 5 层。由此可判断论文的实际组织方式是：

$$
\boxed{
N_{\rm shell}\text{ 沿一条 }M_0\text{ 轨迹固定，}
\quad R_i(t)\text{ 随吸积史变化。}
}
$$

这属于**由表格反推的实现规则**，不是正文显式定义。若改成每个红移都按瞬时 $M(z)$ 重新选择壳层数，就会出现层的生成、合并以及样本/累计量如何映射的问题；论文没有给出这种动态重分层算法。

### 4.2 壳层边界

对于 $N_{\rm shell}$ 个层，论文公式 (27) 使用对数式边界：

$$
R_i(t)=
\left\{
\exp\left[
\frac{i}{N_{\rm shell}}
\ln\left(1+\frac{R_{\rm vir}(t)}{1\,{\rm pc}}\right)
\right]-1
\right\}{\rm pc},
\quad i=0,\ldots,N_{\rm shell}.
\tag{27}
$$

所以 $R_0=0$，$R_{N_{\rm shell}}=R_{\rm vir}$。编号 $i$ 表示同一条质量轨迹中的层级身份，但其物理边界随 $R_{\rm vir}(t)$ 膨胀或收缩。

以 $M_0=10^{12}M_\odot$ 为例，论文表 I 在 $z=0$ 给出约 $R_{\rm vir}=211$ kpc，并列出十层边界：

$$
0, 2.41, 10.61, 38.57, 133.83, 458.46, 1564.65,
5334.13, 18179.03, 61949.37, 211101.44\ {\rm pc}.
$$

### 4.3 壳层代表点与密度

论文不使用壳层体积平均密度，而是在几何中点处评价 NFW 局部密度：

$$
r_{i,\rm mid}(t)=\frac{R_i(t)+R_{i+1}(t)}{2},
$$

$$
\rho_i(t)=
\rho_{\rm NFW}\!\left[r_{i,\rm mid}(t),t\right].
\tag{28}
$$

这是后续硬化率中的局部环境密度。与此同时，给样本重加权时需要的是壳层所含总质量：

$$
M_i(t)=M_{\rm NFW}(<R_{i+1},t)-M_{\rm NFW}(<R_i,t).
$$

后一个公式是**由 NFW 结构和“壳层真实双星数”直接推出的实现量**；论文没有把它单独编号，但如果用中点密度乘壳层体积代替，就不再是严格的 NFW 壳层质量。

### 4.4 局部速度、硬双星边界

论文公式 (24) 定义

$$
v_{{\rm disp},i}(t)
=\sqrt{\frac{2G M(<r_{i,\rm mid},t)}{r_{i,\rm mid}(t)}}.
\tag{24}
$$

当地硬双星边界为

$$
a_{h,i}(t)=\frac{Gm_1}{4v_{{\rm disp},i}^2(t)}.
\tag{23}
$$

因此局部包围质量越大、速度越高，$a_h$ 越小，能够被视为硬双星的初始样本比例也越低。

这里存在一个必须保留的论文内部复现问题：公式 (24) 明写的是 $\sqrt{2GM/r}$，但表 I 与表 III 的若干数值组合更接近 $\sqrt{GM/(2r)}$。后二者速度相差 2 倍，导致 $a_h$ 相差 4 倍。除非能从作者代码或补充材料确认，否则实现应提供明确的 velocity-convention 选项，并把公式分支与表格诊断分支分开，不能静默修改公式。

### 4.5 环境更新时间

halo 密度和速度每

$$
\Delta t=200\ {\rm Myr}
$$

更新一次。在一个全局时间片内部，$R_i,\rho_i,v_{{\rm disp},i},a_{h,i}$ 应视为冻结的背景；轨道则用更小步长演化。论文没有描述在 200 Myr 内对 halo 环境做插值。

论文还假定双星在整个模拟中留在给定壳层，不做壳层间迁移。这里更准确的含义是保留壳层编号；由于 $R_i(t)$ 本身随 halo 演化，双星的代表半径也会随该编号的边界变化。

---

## 5. 双星初始条件：从 $P(a,j)$ 到 $(a_0,e_0)$

### 5.1 物理假设

论文假定 PBH 双星在 halo 外已经形成并保持引力束缚，在进入 halo 前不受 halo 环境影响。进入某个 halo 时间片时，其初态仍服从物质—辐射相等时期得到的早期双星分布。

定义无量纲角动量

$$
j=\sqrt{1-e^2},
\qquad 0<j<1,
$$

以及平均 PBH 间距

$$
\bar x=
\left(\frac{3m}{4\pi f\rho_{\rm eq}}\right)^{1/3}.
$$

附录 D 给出

$$
P(j)=\frac{y^2}{j(1+y^2)^{3/2}},
$$

$$
y=\frac{j}
{0.5\sqrt{1+\sigma_{\rm eq}^2/f^2}\,[x/\bar x]^3},
$$

以及联合分布

$$
P(a,j)=
\frac{3}{4}a^{-1/4}
\left(\frac{f}{\alpha\bar x}\right)^{3/4}
P(j)
\exp\left[-\left(\frac{x(a)}{\bar x}\right)^3\right],
$$

$$
x(a)=
\left(\frac{3am}{4\pi\alpha\rho_{\rm eq}}\right)^{1/4},
\qquad \alpha=0.1.
$$

抽样后通过

$$
e_0=\sqrt{1-j_0^2}
$$

得到偏心率。

### 5.2 抽样范围与论文的描述

论文给出的半长轴抽样区间为

$$
10^{-6}\ {\rm pc}\le a\le1\ {
m pc},
$$

并称使用逆抽样方法。附录文字描述的顺序是先抽 $j$，再由联合分布抽 $a$。但对于真正的联合分布，顺序抽样必须使用相应的边缘分布与条件分布，例如

$$
a\sim P_a(a),
\qquad
j\sim P(j\mid a),
$$

或等价的反向分解。只把联合密度直接当作第二个变量的一维密度，可能改变目标分布。

因此本文档把两件事严格区分：

- **论文原文算法**：按附录 D 的文字顺序复现，并记录其归一化方式；
- **本项目严格联合抽样**：先抽边缘分布，再抽条件分布，用于数学上自洽的基线和对照。

附录 D 的正文称示例抽取 $10^4$ 个样本，而图 19 图注写 $10^5$，两者也不一致。这个数量只是初始分布展示，不要与第 IV 节生产计算中每壳层、每时间片的 $2\times10^6$ 个样本混淆。

### 5.3 硬双星筛选

对每个壳层、每个全局时间片，样本进入轨道积分前必须满足

$$
a_0\le a_{h,i}(t_k).
$$

不满足条件的宽双星不进入该硬双星演化方程。Table II 给出的“硬双星百分比”正是初始分布通过这个局部阈值后的比例。它随 halo 质量、红移和壳层变化，不是一个全局常数。

实现时必须同时保存：

- 总抽样数 $N_{\rm sample}$；
- 通过硬筛选的数量 $N_{\rm hard}$；
- 实际并合数量 $N_{\rm merger}$。

否则无法判断最终的并合比例分母究竟包含全部初始双星，还是只包含硬双星。公式 (32) 使用 $N_{\rm sample}$ 作分母，因此最稳妥的数据结构是让被筛掉的样本自然贡献零次并合，而不是把它们从分母中删除。

---

## 6. 单个硬双星的演化方程

### 6.1 半长轴演化

论文公式 (19) 为

$$
\frac{da}{dt}
=-
\frac{GH\rho_{\rm env}}{v_{\rm disp}}a^2
-\frac{64}{5}
\frac{G^3}{c^5a^3}
(m_1+m_2)(m_1m_2)F(e),
\tag{19}
$$

其中

$$
F(e)=
\frac{1+\frac{73}{24}e^2+\frac{37}{96}e^4}
{(1-e^2)^{7/2}}.
\tag{20}
$$

第一项是平均化的三体硬化，第二项是 Peters 引力波反作用。两个项都使 $a$ 减小，但主导区域不同：

- 较大的硬双星通常先由环境散射收缩；
- $a$ 足够小或 $e$ 足够高后，$a^{-3}F(e)$ 使引力波项迅速接管；
- 最终进入 runaway 式并合。

硬化系数采用

$$
H=14.55\left(1+0.287\frac{a}{a_h}\right)^{-0.95}.
\tag{22}
$$

由于 $H$ 依赖当前 $a/a_h$，即使环境在一个全局时间片中冻结，$H$ 仍应在每个局部积分步重新计算。

### 6.2 偏心率演化

论文公式 (25) 为

$$
\frac{de}{dt}
=\frac{GHK\rho_{\rm env}}{v_{\rm disp}}a
-\frac{304}{15}
\frac{G^3}{c^5a^4}
(m_1+m_2)(m_1m_2)D(e),
\tag{25}
$$

$$
D(e)=
\frac{e+\frac{121}{304}e^3}
{(1-e^2)^{5/2}}.
\tag{26}
$$

第一项描述 binary-single 散射对偏心率的平均改变；第二项描述引力波圆化。论文说 $K$ 取自 Sesana et al. (2006) 的公式 (18)，但目标论文本身没有重写其完整函数和插值实现。因此 $K$ 是代码实现前必须从原始参考文献补齐的外部输入，不能用一个未注明来源的常数代替。

### 6.3 两条竞争时间尺度

理解模拟的最好方式，是比较两个半长轴变化时间尺度：

$$
t_{\rm 3b}\sim
\frac{a}{|\dot a_{\rm 3b}|}
=\frac{v_{\rm disp}}
{GH\rho_{\rm env}a},
$$

$$
t_{\rm GW}\sim
\frac{a}{|\dot a_{\rm GW}|}
\propto
\frac{a^4}{m_1m_2(m_1+m_2)F(e)}.
$$

三体硬化擅长把较宽的硬双星推向更小的 $a$，而引力波在小 $a$、高 $e$ 处变得极强。第 IV 节模拟的物理核心正是判断：**给定 halo 局部环境能否在有限宇宙时间内把一部分早期双星送入引力波主导区。**

附录 E 另给出在 $a=a_h$ 附近的相互作用率估计：

$$
R_{3b}=\frac{2\pi Gm_T n(r,z)a}{v_{\rm disp}},
$$

等质量时 $m_T=3m$、$n=\rho_{\rm env}/m$，所以

$$
R_{3b}(a_h)=
\frac{6\pi G\rho_{\rm NFW}a_h}{v_{\rm disp}},
\qquad
\tau_{3b}=R_{3b}^{-1}.
$$

这个量是环境有效性的诊断，不替代公式 (19)、(25) 的实际积分。

---

## 7. 时间离散与数值积分

### 7.1 两级时间步

论文采用：

$$
\Delta t_{\rm halo}=200\ {\rm Myr},
\qquad
\delta t_{\rm orbit}=2\ {\rm Myr}.
$$

因此每个 halo 时间片通常包含

$$
N_{\rm local}=\frac{200}{2}=100
$$

个局部轨道步。

逻辑上，全局步和局部步职责不同：

- 全局步更新 $M,R_{\rm vir},C,R_i,\rho_i,v_i,a_{h,i}$；
- 局部步在冻结的环境中更新每个样本的 $a,e$。

### 7.2 显式 Euler 更新

论文明确使用 Euler 方法：

$$
a_{n+1}=a_n+delta t\,f_a(a_n,e_n;\rho_i,v_i,a_{h,i}),
$$

$$
e_{n+1}=e_n+delta t\,f_e(a_n,e_n;\rho_i,v_i,a_{h,i}).
$$

其中 $f_a,f_e$ 分别是公式 (19)、(25) 的右端。对一个批次的全部样本，这两个更新应向量化执行。

### 7.3 并合判据

论文只说当 $a$ 接近零时记为并合，没有给出：

- 明确的 $a_{\rm merge}$ 阈值；
- Euler 一步越过 $a=0$ 时如何定位并合时间；
- $e\ge1$、$e<0$ 或浮点溢出时如何处理；
- 是否在 GW 时间尺度短于剩余局部步时直接判定并合；
- 是否使用自适应步长处理临近并合的刚性。

所以“使用 2 Myr Euler 步”可以严格复现，但“稳定而准确地判定并合”仍需要项目自己定义规则。任何钳位、子步进或解析 Peters 剩余时间判据，都应标为数值改进分支，并与原始 Euler 分支对照。

每个样本一旦判为并合，应立即停止继续更新，并且在该批次中只计数一次。

---

## 8. 一个时间片内究竟做什么

对给定 $M_0$、时间片 $t_k\to t_{k+1}$ 和壳层 $i$，完整操作是：

1. 由 $M(z_k)$ 和浓度模型构造当前 NFW halo；
2. 由公式 (27) 更新该层边界和代表半径；
3. 求 $\rho_i,v_i,a_{h,i}$；
4. 从 $P(a,j)$ 抽取 $N_{\rm sample}$ 组初态；
5. 转换为 $e_0=\sqrt{1-j_0^2}$；
6. 用 $a_0\le a_{h,i}$ 标记硬双星；
7. 对硬双星做最多 100 个 Euler 局部步；
8. 记录在这 200 Myr 内并合的样本数 $N_{{\rm merger},i,k}$；
9. 由 NFW 壳层质量求真实双星数 $N_{{\rm BBH},i,k}$；
10. 把样本并合比例重加权成该层的物理并合率。

论文生产计算称每个壳层、每个时间片都取

$$
N_{\rm sample}=2\times10^6.
$$

这会让内外壳层拥有近似相同的抽样噪声预算。未经重加权时，内层通常因密度高而产生更多并合；重加权后，外层可能因为体积和真实双星数更大而占据重要甚至主导贡献。

### 8.1 样本是否跨全局时间片继承

这是论文算法中最重要的未决点之一。正文说每个时间步选择一批初始 $(a_0,e_0)$，并在 $t$ 到 $t+\Delta t$ 之间演化；又说每个壳层、每个时间步使用 $2\times10^6$ 个样本。这更接近“每个环境时间片独立抽取一批早期双星并估计该时间片的并合概率”。

但论文没有说明：

- 上一时间片未并合的 Monte Carlo 个体是否传入下一时间片；
- 新吸积进入 halo 的双星数如何与原有存量区分；
- $N_{{\rm BBH},i,k}$ 是否扣除过去已并合的系统；
- 壳层随增长容纳的新质量如何分配给 cohort。

因此，严格复现时应把“独立时间片重新抽样”和“cohort 跨时间继承”作为两个不同算法。就现有文字而言，前者更贴近论文公开描述，但不能把它声称为作者代码已经确认的事实。

---

## 9. 从样本并合数到真实 halo 并合率

### 9.1 真实双星数

设壳层 NFW 质量为 $M_i(t_k)$。若 $f_{\rm PBH}=1$，并且 halo 质量中有比例 $f_{\rm binary}$ 位于 PBH 双星内，则由质量守恒直接得到

$$
N_{{\rm BBH},i,k}
=\frac{f_{\rm binary}M_i(t_k)}{2m}.
$$

一般化后为

$$
N_{{\rm BBH},i,k}
=\frac{f_{\rm PBH}f_{\rm binary}M_i(t_k)}{2m}.
$$

这是**由论文比例设定直接推出的系统计数公式**；论文第 IV 节没有明确写出分母到底采用 $m$ 还是 $2m$。实现时必须明确 $f_{\rm binary}$ 是“双星中所含 PBH 的质量比例”还是“PBH 对象比例”，否则会产生因子 2 的歧义。

### 9.2 Monte Carlo 重加权

样本在一个时间片中的并合概率估计为

$$
\hat p_{i,k}
=\frac{N_{{\rm merger},i,k}}{N_{\rm sample}}.
$$

该壳层在该时间片的期望并合数为

$$
\widehat{\Delta N}_{{\rm merge},i,k}
=N_{{\rm BBH},i,k}\hat p_{i,k},
$$

相应的平均率为

$$
\widehat R_{i,k}
=\frac{N_{{\rm BBH},i,k}}{N_{\rm sample}}
\frac{N_{{\rm merger},i,k}}{\Delta t}.
$$

对所有壳层相加得到逐时间片的单 halo 率：

$$
\boxed{
R_{\rm halo}(M_0,t_k)
=\sum_i
\frac{N_{{\rm BBH},i,k}}{N_{\rm sample}}
\frac{N_{{\rm merger},i,k}}{\Delta t}
}
$$

其单位可转换为 ${\rm yr}^{-1}\,{\rm halo}^{-1}$。

### 9.3 公式 (32) 的时间求和歧义

论文公式 (32) 写成对 $i,t$ 都求和的形式：

$$
R_{\rm halo}
=\sum_i\sum_t
\frac{N_{{\rm PBH\ binaries},i,t}}{N_{\rm sample}}
\frac{N_{{\rm merger},i,t}}{\Delta t}.
\tag{32}
$$

若要画“随红移变化的并合率”，不能在每个红移点再把全部 $t$ 求和；应使用上面的逐时间片形式。若对 $k$ 累加，则得到的是所有时间片贡献的总和，或在乘回 $\Delta t$ 后得到累计并合数：

$$
N_{\rm merge,cum}
=\sum_k R_{\rm halo}(t_k)\Delta t.
$$

所以代码应分别保存：

- `rate_per_time_bin`：逐红移/逐时间片率；
- `expected_mergers_per_bin`：该片的期望并合数；
- `cumulative_mergers`：前两者积分得到的累计量。

不能让一个名为 `R_halo` 的标量同时承担三种意义。

### 9.4 Monte Carlo 误差

对于独立 Bernoulli 计数，未加权的并合比例误差约为

$$
\sigma_{p,i,k}
\simeq
\sqrt{\frac{\hat p_{i,k}(1-\hat p_{i,k})}{N_{\rm sample}}}.
$$

当外层只有极少样本并合时，重加权后会出现明显尖峰。论文对此做了多项式平均，并丢弃第一个时间点，因此部分图从约 $z\simeq10$ 开始，而非模拟起点 $z=12$。但论文没有公开多项式阶数、窗口、权重或边界处理，所以平滑后的曲线不能仅凭文字逐点复现。代码应优先保存原始计数和未平滑率，平滑只作为派生输出。

---

## 10. 从单 halo 率到宇宙共动并合率

### 10.1 halo 质量网格

论文使用 50 个今天质量点，按对数覆盖

$$
10^3M_\odot\le M_0\le10^{15}M_\odot.
$$

对每条质量轨迹都执行前述壳层—时间—Monte Carlo 计算，得到 $R_{\rm halo}(M,z)$。

### 10.2 HMF 积分

共动体积并合率的基本形式是

$$
\boxed{
\mathcal R(z)
=\int dM\,
\frac{dn(M,z)}{dM}
R_{\rm halo}(M,z)
}
$$

单位为

$$
{\rm Gpc}^{-3}\,{\rm yr}^{-1}.
$$

论文基准采用 Press–Schechter halo mass function，并用 Jenkins HMF 做敏感性比较。论文报告在 $z<1$ 时 Press–Schechter 给出的结果可高约 6–8 倍，较高红移处差距缩小。这说明 HMF 不是外围输入，而是总体率归一化的重要系统误差源。

实现时还需明确积分变量是目标红移处的瞬时 halo 质量 $M(z)$，还是由今天质量 $M_0$ 标记的轨迹。论文用 $M_0$ 组织质量吸积轨迹，但没有详述把轨迹网格插入目标红移 HMF 的全部数值细节；这需要在接口中显式保留 `track_mass_M0` 与 `instantaneous_mass_Mz` 两个字段。

### 10.3 与直接俘获通道相加

第 IV 节 binary-single 结果最终与前文的直接双体俘获率相加：

$$
\mathcal R_{\rm total}
=\mathcal R_{\rm capture}
+\mathcal R_{\rm binary-single}.
$$

这只是两个通道率的直接相加。论文没有建立共享 PBH 库，也没有处理一个 PBH 已在某通道消耗后不能再进入另一个通道的联合耗尽问题。因此它是“稀疏事件、背景固定”近似下的总率，而不是两通道完全自洽的总体演化。

论文还给出总率随 $f_{\rm PBH}^2$ 缩放的讨论。该缩放是在其固定环境与比例设定下的重标度，不等价于重新计算不同 $f_{\rm PBH}$ 对 halo 结构、双星初始分布和耗尽历史的反馈。

---

## 11. 从输入到输出的完整伪代码

下面的伪代码严格区分“轨迹标签”“瞬时 halo 状态”“Monte Carlo 样本”和“物理重加权量”：

```text
INPUT:
    cosmology
    present-day halo masses M0_grid (50 log-spaced points, 1e3--1e15 Msun)
    mass-accretion model M(z | M0)
    concentration model C(M,z)             # paper: Ludlow16
    HMF dn/dM                              # paper baseline: Press-Schechter
    PBH masses m1 = m2 = 30 Msun
    f_PBH = 1, f_binary = f_single = 0.5
    initial binary sampler P(a,j)
    dt_global = 200 Myr
    dt_local  = 2 Myr
    N_sample  = 2e6 per shell per global time bin

FOR each M0 in M0_grid:
    choose start redshift
        if M0 >= 5e4 Msun: z_start = 12
        else: lower z_start until halo contains about 30 PBHs

    choose N_shell from the M0 track category
    build cosmic-time bins from z_start to z=0

    FOR each global bin k = [t_k, t_k + dt_global]:
        z_k = z(t_k)
        M_k = M(z_k | M0)
        C_k = concentration(M_k, z_k)
        construct NFW(M_k, C_k, z_k)
        compute Rvir_k
        compute shell boundaries R_i(k) from Eq. (27)

        FOR each shell i:
            r_mid = (R_i + R_{i+1}) / 2
            rho_i = rho_NFW(r_mid)
            Menc_i = M_NFW(<r_mid)
            v_i = sqrt(2 G Menc_i / r_mid)       # printed Eq. (24)
            ah_i = G m1 / (4 v_i^2)              # Eq. (23)

            shell_mass = M_NFW(<R_{i+1}) - M_NFW(<R_i)
            N_BBH_actual = f_PBH * f_binary * shell_mass / (2m)
                # factor convention must be documented

            sample N_sample independent (a0,j0) from P(a,j)
            e0 = sqrt(1-j0^2)
            hard_mask = (a0 <= ah_i)

            a = a0[hard_mask]
            e = e0[hard_mask]
            merged = false for every retained sample

            REPEAT up to dt_global/dt_local = 100 times:
                H = H(a/ah_i)
                K = Sesana-fit K(...)
                da_dt = three_body_hardening + GW_shrinkage
                de_dt = three_body_eccentricity + GW_circularization
                a <- a + dt_local * da_dt
                e <- e + dt_local * de_dt
                identify newly merged samples
                stop updating merged samples

            N_merger[i,k] = count(merged)
            p_merge[i,k] = N_merger[i,k] / N_sample
            expected_mergers[i,k] = N_BBH_actual * p_merge[i,k]
            rate_shell[i,k] = expected_mergers[i,k] / dt_global

        rate_halo[M0,k] = sum_i rate_shell[i,k]
        cumulative_halo[M0,k] = time integral of rate_halo

AT each requested redshift z_k:
    interpolate rate_halo onto the HMF mass variable consistently
    rate_comoving[z_k] = integral dM (dn/dM)(M,z_k) * rate_halo(M,z_k)

OUTPUT:
    halo histories and shell diagnostics
    raw sample counts N_hard, N_merger
    shell rate, per-halo rate, cumulative mergers
    comoving binary-single merger rate
    optional total = binary-single + direct-capture rate
```

---

## 12. 每一级输出是什么，它回答什么问题

### 12.1 halo/壳层诊断输出

建议最先保存：

- $M(z),R_{\rm vir}(z),C(z)$；
- 每层 $R_i,r_{\rm mid},M_i,\rho_i,v_i,a_{h,i}$；
- $N_{\rm shell}$ 及其选择依据；
- 所用速度约定和浓度模型名称。

这些量回答：“双星被放进了怎样的环境？”论文图 10、11 主要属于这一层。

### 12.2 样本级输出

- 初始 $a_0,j_0,e_0$；
- `hard_mask`；
- 并合标记与并合所在局部步；
- 必要时保存少量轨道的 $a(t),e(t)$，而不是保存全部 $2\times10^6$ 条轨迹。

这些量回答：“Monte Carlo 样本为何被保留，以及通过哪种动力学走向并合？”

### 12.3 统计输出

- $N_{\rm sample},N_{\rm hard},N_{\rm merger}$；
- 原始 $hat p_{i,k}$ 与统计误差；
- $N_{{\rm BBH},i,k}$；
- 壳层重加权前、后的贡献；
- 未平滑率和平滑率。

这些量回答：“图中的并合率究竟来自高并合概率，还是来自巨大的真实双星数量？”论文图 12–14 主要展示这一层。

### 12.4 最终物理输出

- $R_{\rm halo}(M_0,z)$，单位 ${\rm yr}^{-1}\,{\rm halo}^{-1}$；
- $\mathcal R_{\rm binary-single}(z)$，单位 ${\rm Gpc}^{-3}\,{\rm yr}^{-1}$；
- 不同 HMF、浓度关系和 $f_{\rm PBH}$ 的分支；
- 与直接俘获通道相加后的 $\mathcal R_{\rm total}(z)$。

论文图 15–17 主要属于这一层。

---

## 13. 最容易混淆的逻辑关系

### 13.1 壳层编号固定，不等于壳层物理结构固定

沿 $M_0$ 轨迹固定的是层数和层级编号；随吸积史变化的是 $M(z),R_i(t),\rho_i(t),v_i(t),a_h(t)$。如果把所有壳层量冻结在 $z=0$，就失去了论文研究 halo 演化的核心。

### 13.2 局部密度不等于壳层总质量

$\rho_{\rm NFW}(r_{\rm mid})$ 进入硬化微分方程；壳层积分质量 $M_i$ 决定真实双星数。两者职责不同，不能共用一个近似量而不说明。

### 13.3 硬双星比例不等于并合概率

$a_0<a_h$ 只说明样本可用硬化公式继续演化；能否在 200 Myr 内并合还取决于 $a_0,e_0,\rho,v,H,K$ 以及 GW 项。

### 13.4 原始并合数不等于物理壳层贡献

相同 $N_{\rm sample}$ 下，内层原始并合数常较高；但外层的 $N_{\rm BBH}$ 可以大得多。只有完成公式 (32) 的人口重加权，才能比较壳层的物理贡献。

### 13.5 单 halo 率不等于宇宙总体率

高质量 halo 单体率可能高，但数量稀少；低质量 halo 单体率低，却可能数量庞大。HMF 积分决定最后哪个质量区间主导宇宙总体率。

---

## 14. 论文没有完全封闭的实现细节

在开始正式编码前，应把下列项目作为显式配置或待核对项：

1. **速度约定冲突**：公式 (24) 与表 I/III 数值并不完全一致；
2. **$K$ 的实现**：需从 Sesana et al. (2006) 补齐函数、参数范围和插值；
3. **联合分布抽样**：附录 D 的文字顺序是否真正等价于目标 $P(a,j)$；
4. **生产样本的继承关系**：每个 200 Myr 独立重采样，还是继承未并合 cohort；
5. **并合阈值**：何谓 $a\to0$，以及 Euler 越界如何处理；
6. **偏心率边界**：数值上如何保证 $0\le e<1$；
7. **实际双星数定义**：$f_{\rm binary}$ 的质量比例/对象比例及因子 2；
8. **公式 (32) 的时间求和**：逐时间片率与累计量必须分开；
9. **平滑细节**：多项式阶数、拟合窗口和权重未公开；
10. **HMF 的质量变量**：瞬时 $M(z)$ 与轨迹标签 $M_0$ 的插值映射；
11. **逃逸与耗尽**：论文忽略，不能在解释结果时当作已经自洽处理；
12. **最后一个不足 200 Myr 的时间片**：论文未说明舍弃、缩短还是延伸。

这些并不妨碍建立主程序骨架，但每一项都应在结果元数据中留下选择记录，否则不同实现可能产生“都声称按论文、实际算法却不同”的情况。

---

## 15. 与本项目后续接口的对应关系

为了让后续程序既能复现论文，也能替换 Prada12-HMF 等模型，建议保持以下责任边界：

```text
mass_history(M0, z, model)
    -> instantaneous M(z)

concentration(M, z, model)
    -> C(M,z)

halo_structure(M0, z, history_model, concentration_model,
               shell_policy, velocity_convention)
    -> Rvir, NFW parameters, shell boundaries,
       shell mass, rho_mid, v_disp, a_h

initial_binary_sampler(N, sampler_mode)
    -> a0, j0, e0, sampling weights

evolve_binary_batch(a0, e0, environment, dt_global, dt_local,
                    K_model, merger_criterion)
    -> hard mask, merger mask, merger time, diagnostics

reweight_shell_statistics(counts, actual_binary_population, dt_global)
    -> shell rate and uncertainty

integrate_halo_rate(shell_rates)
    -> per-halo rate versus time/redshift

integrate_hmf(per_halo_rate, dn_dM)
    -> comoving merger-rate density
```

其中：

- `history_model="paper"` / `concentration_model="Ludlow16"` 对应论文分支；
- `history_model="Prada12-HMF"` 对应本项目后续默认分支；
- `velocity_convention="eq24"` 应严格实现印刷公式；
- 表格匹配只能作为单独诊断分支；
- `sampler_mode="paper_order"` 与 `sampler_mode="strict_joint"` 必须分开；
- 原始计数永远先保存，HMF 积分和平滑不得覆盖它们。

---

## 16. 最终理解：这个模拟真正做了什么

从物理上说，它研究的是一条“环境催化”链：

$$
\text{halo 吸积增长}
\rightarrow
\text{局部密度和速度改变}
\rightarrow
\text{硬双星阈值与三体硬化效率改变}
\rightarrow
\text{双星进入 GW 主导区的概率改变}
\rightarrow
\text{单 halo 并合率改变}
\rightarrow
\text{宇宙总体并合率改变}.
$$

从数值上说，它做的不是“真实生成一个 halo 中的全部 PBH 并逐个演化”，而是：

$$
\text{在规则化环境网格上估计条件并合概率，}
$$

再用实际人口和 HMF 两次重加权：

$$
\underbrace{\hat p_{\rm merge}(M,z,i)}_{
\text{轨道 Monte Carlo}}
\times
\underbrace{N_{\rm BBH}(M,z,i)}_{
\text{壳层人口重加权}}
\times
\underbrace{\frac{dn}{dM}(M,z)}_{
\text{halo 丰度重加权}}.
$$

这三个因子分别回答：

1. 一个双星在该环境中多容易并合？
2. 这种环境中实际有多少双星？
3. 宇宙中实际有多少这样的 halo？

只有三者都定义清楚，最后的 ${\rm Gpc}^{-3}\,{\rm yr}^{-1}$ 才具有明确意义。这也是后续编写程序时最应该保持的数据流主线。
