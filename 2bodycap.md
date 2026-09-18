# 两体直接捕获总并合率与 Fig. 8 排查记录

## 0. 先概括：问题究竟出在哪里

原来的 Fig. 8 在低红移附近看起来尚可，但到 \(z=12\) 时只有约
\(10\!-\!15\ {\rm Gpc}^{-3}{\rm yr}^{-1}\)，而论文原图约为
\(130\!-\!160\ {\rm Gpc}^{-3}{\rm yr}^{-1}\)，相差约一个数量级。

排查后的结论是：**主要问题不在 HMF 的质量上下界，也不在 \(h\) 的单位换算、
\(dM\) 与 \(d\ln M\) 的 Jacobian、两体捕获截面或 NFW 径向积分；真正造成高红移
差异的，是“严格的 Eq. (18) 积分”与“论文 Fig. 8 实际采用的数值轨迹/网格约定”
并不相同。** 此外，Prada12 分支还把 Eq. (23) 的闭式 \(\sigma(M)\) 近似用于了
论文明确引用 Eqs. (12)--(22) 的场合，导致低质量晕浓度和总率出现额外偏差。

本次找到的具体问题依次为：

1. 旧 Fig. 8 直接在每个红移的当前晕质量网格上计算严格 Eq. (18)。这个积分在
   数学上自洽，但不能复现论文给出的高红移曲线。
2. 论文原始矢量曲线与以下算法高度一致：从今天的 \(M_0\) 出发演化 50 条质量
   吸积轨迹，再把非均匀的轨迹质量与另一条均匀对数 HMF 网格按数组下标配对，
   最后对 \(\ln M_0\) 积分。
3. 这种“按数组下标配对”不是严格 Eq. (18)，更像作者实现中的数值网格约定。
   论文没有公开 Fig. 8 的源代码，因此这是由原始矢量数据支持的反推结论，不能
   当成新的物理效应。
4. Prada12 应使用 WMAP5 线性功率谱计算 \(\sigma(M,z)\)，再代入 Eqs. (12)--(22)；
   旧实现使用 Eq. (23) 的闭式近似，在低质量端产生了明显偏差。
5. Ludlow16 的公开 Appendix C 质量吸积史与论文随附的 `MAH_corea.pdf` 矢量轨迹
   并非逐点一致。为了复现 Fig. 8，Ludlow 分支采用该中间图反求的
   \(\alpha(M_0)\)、\(\beta(M_0)\)，但严格的公开公式仍保留在原模块中。

因此现在代码明确保留两套入口：

- `comoving_capture_rate_eq18_gpc3_per_year`：严格的当前质量 Eq. (18) 积分；
- `comoving_capture_rate_with_history_gpc3_per_year`：接收外部质量吸积史函数，
  按 Fig. 8 的轨迹和网格约定完成总率积分。

二者不能混为一谈，也不能用 Fig. 8 的实现替代一般物理计算。

---

## 1. 被排查的完整计算链

两体直接捕获总率的程序调用链是

\[
\text{PBH 质量分布}
\longrightarrow \langle\sigma_{\rm cap}v\rangle
\longrightarrow R_{\rm halo}(M,z)
\longrightarrow \frac{dn}{d\ln M}(M,z)
\longrightarrow \mathcal R(z).
\]

更具体地说：

1. `pbhmassfunction.py` 给出离散或连续的 PBH 质量概率权重；
2. `twobodycapture.py` 计算引力波耗散两体捕获截面和质量双重积分；
3. `haloconcentration.py` 给出 Ludlow16 或 Prada12 浓度；
4. `halostructure.py` 由 \((M,z,C)\) 构造 \(R_{\rm vir}\)、\(R_s\)、\(\rho_s\)、
   \(\rho_{\rm NFW}(r)\) 和截断速度分布；
5. `capture_rate_per_halo_a14_per_year` 实现 Appendix A14 的每晕捕获率；
6. `press_schechter_halo_mass_function_hmfcalc` 给出 \(dn/d\ln M\)；
7. 最后对晕质量积分得到单位为
   \({\rm Gpc}^{-3}{\rm yr}^{-1}\) 的共动总率。

这条链很长，所以排查时不能看到总率不对就直接调整归一化。必须先确认究竟是
单晕率、HMF、质量轨迹还是最终积分测度出了问题。

---

## 2. 第一项检查：HMF 的质量上下界是否设置错误

### 2.1 `hmf` 原生质量约定

本工程公开接口使用物理质量 \(M_{\rm phys}\)，单位为 \(M_\odot\)。Python
`hmf.MassFunction` 的内部约定则是

\[
M_{\rm min}^{\rm hmf}
=\log_{10}\!\left(\frac{M}{M_\odot/h}\right),
\qquad
M_{\rm max}^{\rm hmf}
=\log_{10}\!\left(\frac{M}{M_\odot/h}\right).
\]

物理质量与 `hmf.m` 数值之间满足

\[
M_{\rm hmf}=hM_{\rm phys},
\qquad
M_{\rm phys}=\frac{M_{\rm hmf}}{h}.
\]

因此，当公开接口要求

\[
10^3M_\odot\le M_{\rm phys}\le10^{15}M_\odot
\]

时，传给 `MassFunction` 的正确参数是

\[
M_{\min}=\log_{10}(10^3h),
\qquad
M_{\max}=\log_{10}(10^{15}h),
\]

而不是直接把 `1e3`、`1e15` 传给 `Mmin`、`Mmax`。

当前函数已经正确完成了这一步换算：

```python
Mmin=np.log10(minimum_halo_mass_msun * hubble_h)
Mmax=np.log10(maximum_halo_mass_msun * hubble_h) + 0.5 * dlog10m
```

这里给 `Mmax` 增加半个步长，是因为 `hmf` 的质量数组由类似
`np.arange(Mmin, Mmax, dlog10m)` 的右开区间产生。函数随后再用物理质量范围
筛选，因此不会把额外点留在返回结果中。

### 2.2 数密度单位换算

`MassFunction.dndlnm` 的原始单位为

\[
h^3\,{\rm Mpc}^{-3}.
\]

公开接口返回物理 \({\rm Mpc}^{-3}\)，所以代码执行

\[
\left.\frac{dn}{d\ln M}\right|_{\rm physical}
=h^3
\left.\frac{dn}{d\ln M}\right|_{\rm hmf}.
\]

由于

\[
d\ln(M/h)=d\ln M,
\]

质量只差一个常数因子时不产生新的 Jacobian。

### 2.3 实际检查结果

使用

```python
minimum_halo_mass_msun = 1.0e3
maximum_halo_mass_msun = 1.0e15
dlog10m = 0.1
```

实际得到：

- 网格点数：121；
- 第一项：\(999.9999999999999M_\odot\)；
- 最后一项：\(1.0000000000000244\times10^{15}M_\odot\)；
- 所有返回的 \(dn/d\ln M\) 均为有限正数；
- `dndlnm = m*dndm` 的最大相对误差为 0；
- 用 \(d\ln M\) 和 \(d\log_{10}M\) 两种等价形式积分，结果一致。

结论：**`minimum_halo_mass_msun=1e3` 与
`maximum_halo_mass_msun=1e15` 对当前封装函数是正确的，HMF 上下界和
\(h\) 单位不是 Fig. 8 相差一个数量级的原因。**

需要保留一个来源边界：目标论文明确提到 HMFcalc，而本工程调用的是 Python
`hmf` 包来实现同类计算。HMFcalc 当前也以 `hmf` 为计算后端，但作者当时使用的
版本、宇宙学参数、转移函数和全部界面选项没有公开，所以不能声称两边配置完全
相同。

---

## 3. 第二项检查：用 Bird 论文 Fig. 2 验证 HMF 与单晕率链

### 3.1 为什么选择 Bird Fig. 2

`Did LIGO Detect Dark Matter?` 的 Fig. 2 正好画出

\[
\frac{d\mathcal V}{d\ln M}
=R_{\rm halo}(M,0)\frac{dn}{d\ln M}(M,0),
\]

并比较 Tinker、Press--Schechter、Jenkins 等 HMF。它可以同时检查：

- 每晕捕获率的数量级；
- HMF 质量坐标和密度单位；
- \(d\ln M\) 积分测度；
- 低质量截止对总率的影响。

### 3.2 本次复现结果

以 Bird 使用的 \(30M_\odot\) 单色 PBH 和
\(M_{\min}=400M_\odot\) 为例，得到

| HMF | 本次积分结果 / \({\rm Gpc}^{-3}{\rm yr}^{-1}\) |
|---|---:|
| Tinker08 | 2.5406 |
| Press--Schechter | 1.4356 |
| Jenkins | 0.003448 |

Bird 的近似式在该截止处给出 Press--Schechter 约 2，文中给出的 Jenkins 结果约
0.02。数值没有完全相同，尤其 Jenkins 外推仍受具体 HMF 配置影响，但以下结构均
成功复现：

1. 小质量晕主导总率；
2. Tinker 给出最高的低质量贡献；
3. Press--Schechter 比 Tinker 低；
4. Jenkins 在小质量端被压低多个数量级；
5. 改变 \(M_{\min}\) 会强烈改变积分总率。

对应诊断文件为：

- `notebooks/did_ligo_fig2_hmf_diagnostic.ipynb`；
- `did_ligo_fig2_reproduction.png`；
- `did_ligo_fig2_reproduction.csv`。

这个交叉检验说明：HMF 和每晕率的基本量纲、趋势及积分测度没有足以解释 Fig. 8
高红移约 10 倍差异的错误。Bird 复现仍不是逐像素一致，因此它只能用于排除大类
错误，不能反过来证明目标论文的全部 HMFcalc 参数已经确定。

---

## 4. 第三项检查：两体捕获和每晕率是否错误

对目标论文的 Fig. 3 和 Fig. 7 进行过以下比较：

- Ludlow16 的 \(z=0\) 每晕率通常与原始矢量曲线相差约 2%--6%；
- 大多数质量轨迹上的红移演化只相差几个百分点；
- Ludlow16 在 \(z\simeq5.9\) 附近可出现约 12% 的偏差；
- 最高质量轨迹还有已经记录的质量吸积史差异；
- 旧 Prada12 在 \(10^3M_\odot\) 附近可高约 42%，进一步追踪后发现主要来自
  \(\sigma(M,z)\) 近似，而不是捕获截面。

因此，下列部分不是 Fig. 8 高红移低一个数量级的主要原因：

- Appendix A14 的 PBH 双质量积分；
- 唯一配对因子 \(1/2\)；
- \(\langle\sigma_{\rm cap}v\rangle\)；
- NFW 的 \(\rho^2\) 径向积分；
- Bird 截断 Maxwell 三维速度密度的积分测度；
- \({\rm Mpc}^{-3}\to{\rm Gpc}^{-3}\) 的 \(10^9\) 换算。

这一步非常重要：如果每晕率已经在几个百分点内复现，就不应该为了总率曲线而
改动捕获截面或任意乘一个红移相关因子。

---

## 5. Prada12 的具体问题：错误使用 Eq. (23) 近似

### 5.1 两种 Prada12 计算口径

Prada12 的核心关系可以写成

\[
x=\left(\frac{\Omega_{\Lambda0}}{\Omega_{m0}}\right)^{1/3}a,
\qquad
\sigma'=B_1(x)\sigma(M,z),
\]

\[
C(M,z)=B_0(x)\,\mathcal C(\sigma').
\]

其中 \(B_0\)、\(B_1\)、\(\mathcal C\) 对应 Prada12 Eqs. (14)--(22)。这里仍然
需要外部给出线性方差 \(\sigma(M,z)\)。有两种做法：

1. 由线性功率谱、top-hat 窗函数和增长因子真正计算 \(\sigma(M,z)\)；
2. 使用 Prada12 Eq. (23) 给出的闭式近似。

Eq. (23) 适合快速估算和学习公式，但目标论文正文明确说 Prada 曲线采用
Eqs. (12)--(22)，而不是 Eq. (23)。旧 Fig. 8 把闭式近似当成了作者实际使用的
\(\sigma\)，这在低质量端会放大浓度并改变每晕率。

### 5.2 修正方法

`haloconcentration.py` 中新增函数

```python
hc.linear_sigma_hmfcalc_wmap5(mass_msun, z)
```

使用 WMAP5 宇宙学、CAMB 转移函数和 top-hat 滤波计算 \(\sigma(M,z)\)。随后

```python
hc.concentration_prada12_hmf_sigma(mass_msun, z)
```

把该 \(\sigma\) 代入 Prada12 Eqs. (12)--(22)。在高峰高区域，根据目标论文的
说明，把 Prada 上翘分支限制在当时的最小浓度，而不是无限沿 U 形右支增长。

原来的 Eq. (23) 实现没有删除，仍保留在 `haloconcentration.py` 中用于学习和
快速估算。新旧函数代表不同近似层次，不能相互冒充。

---

## 6. Ludlow16 质量吸积史的具体问题

公开 Appendix C 使用

\[
M(z)=M_0(1+z)^\alpha e^{\beta z}.
\]

旧程序完全根据当前浓度模型计算 \(z_{-2}\)、\(\alpha\)、\(\beta\)。这个实现遵守
公开公式，但与目标论文随附 `MAH_corea.pdf` 中的 13 条矢量轨迹存在质量相关的
差异，尤其在高质量和高红移处更明显。

为了复现 Fig. 8，本次从该矢量中反求出以下轨迹参数：

| \(\log_{10}(M_0/M_\odot)\) | \(\alpha\) | \(\beta\) |
|---:|---:|---:|
| 3 | 0.100850 | -0.297334 |
| 4 | 0.108239 | -0.316447 |
| 5 | 0.116726 | -0.339144 |
| 6 | 0.126482 | -0.366252 |
| 7 | 0.137709 | -0.398855 |
| 8 | 0.150990 | -0.439511 |
| 9 | 0.166727 | -0.490914 |
| 10 | 0.185513 | -0.557620 |
| 11 | 0.207982 | -0.647224 |
| 12 | 0.233766 | -0.770728 |
| 13 | 0.270119 | -0.903985 |
| 14 | 0.331111 | -1.013984 |
| 15 | 0.392102 | -1.123982 |

这些数字只用于 `notebooks/2bodycaptot.ipynb` 中 `_figure8_mass_history` 的
Ludlow Fig. 8 分支。它们是从论文提供的中间矢量结果恢复的数据，不是用 Fig. 8
的总率曲线拟合出来的自由参数。

Prada 分支则先用修正后的 Prada 浓度计算 \(C(M_0,0)\)，再代入 Appendix C。
质量吸积史的背景常数采用论文 Appendix C 所属的 Ludlow/Planck 口径和
\(A=798\)。这是目前与原始曲线最一致、且能由论文结构支持的组合，但作者没有
公开完整代码，所以仍应保留“数值复现推断”的标签。

---

## 7. 决定性问题：严格 Eq. (18) 与论文 Fig. 8 网格约定不同

### 7.1 严格的 Eq. (18)

严格积分应在同一个实际当前质量 \(M\) 上同时评价两个因子：

\[
\mathcal R(z)
=\int_{M_{\min}}^{M_{\max}}
R_{\rm halo}(M,z)\frac{dn}{dM}(M,z)\,dM,
\]

等价地，

\[
\mathcal R(z)
=\int
R_{\rm halo}(M,z)\frac{dn}{d\ln M}(M,z)\,d\ln M.
\]

如果使用今天的质量标签 \(M_0\)，变量替换后还必须带上正确的映射和 Jacobian，
不能只把 \(M(z;M_0)\) 塞进其中一部分。

`comoving_capture_rate_eq18_gpc3_per_year` 保留了这种严格的当前质量积分。

### 7.2 从论文矢量曲线反推的作图算法

原始 Fig. 8 可以由以下过程复现到约 10%：

1. 在今天的质量范围
   \(M_0=10^3\!-\!10^{15}M_\odot\) 上建立 50 个对数等距点；
2. 对每个 \(M_0\) 沿质量吸积史计算 \(M(z;M_0)\)；
3. 在这些非均匀的轨迹质量上计算每晕率
   \(R_{\rm halo}[M(z;M_0),z]\)；
4. 另行在
   \([\min M(z),\max M(z)]\) 之间建立 50 个均匀对数 HMF 质量点；
5. 按数组下标把第 \(i\) 个每晕率与第 \(i\) 个 HMF 数密度相乘；
6. 对今天的 \(\ln M_0\) 进行梯形积分。

关键在于：步骤 3 的 `track_masses[i]` 一般不等于步骤 4 的
`paired_hmf_masses[i]`。因此，被相乘的每晕率和 HMF 并不是在同一个实际质量上
求值。

该处理会显著提高高红移结果，因为质量轨迹在高红移下变得非均匀，数组下标配对
改变了低质量端高数密度部分的权重。\(z=0\) 时
\(M(z;M_0)=M_0\)，两条网格重合，所以这个问题在低红移几乎看不出来；随着红移
升高，差异迅速放大。这正好解释了旧图为何在 \(z=0\) 附近正常、在 \(z=12\)
低一个数量级。

### 7.3 分层实验如何锁定该问题

以单色 PBH 为例：

| 算法 | Ludlow16 \(R(12)\) | Prada12 \(R(12)\) |
|---|---:|---:|
| 旧的严格当前质量积分 | 10.38 | 13.35 |
| 追踪 \(M_0\)，但在实际 \(M(z)\) 上严格配对 HMF | 约 71.5 | 约 60 |
| 使用论文式“轨迹网格 + 均匀 HMF 网格按下标配对” | 134.47 | 135.68 |
| 论文原始矢量曲线 | 139.59 | 130.98 |

这个分层实验说明：

- 仅仅加入质量吸积史只能解释一部分差异；
- 真正把结果从约 60--70 推到约 130--140 的，是两套质量网格的下标配对；
- 不需要额外乘一个人为的 10 倍归一化常数。

必须再次强调：**这里复现的是作者画图时很可能采用的数值约定，并不是证明严格
Eq. (18) 应这样计算。** 在后续物理预测中，应优先使用严格接口；只有比较论文
Fig. 8 时才使用论文专用接口。

---

## 8. PBH 质量分布为什么不是问题来源

目标论文 Fig. 8 的对数正态曲线与单色曲线之比，在所有红移处几乎恒定为

\[
\frac{\mathcal R_{\rm lognormal}}
{\mathcal R_{\rm mono}}\simeq1.1352.
\]

当前 Appendix A14 双质量积分也复现了相同的恒定比例。这说明：

- 单色质量函数的调用正确；
- 连续质量分布权重的归一化正确；
- 双质量求积和唯一配对因子没有引入红移相关错误；
- Fig. 8 的主要差异来自晕质量和 HMF 部分，而不是 PBH 质量函数。

---

## 9. 最终修正结果

### 9.1 \(z=12\) 的旧值、原文值和修正值

| 模型与质量函数 | 旧结果 | 论文原图 | 修正结果 |
|---|---:|---:|---:|
| Ludlow16 + 单色 | 10.38 | 139.59 | 134.47 |
| Ludlow16 + 对数正态 | 11.78 | 158.46 | 152.64 |
| Prada12 + 单色 | 13.35 | 130.98 | 135.68 |
| Prada12 + 对数正态 | 15.15 | 148.68 | 154.02 |

单位均为 \({\rm Gpc}^{-3}{\rm yr}^{-1}\)。

### 9.2 全红移范围误差

与从论文原始 PDF 矢量路径恢复的 13 个
\(z=0,1,\ldots,12\) 数据点比较：

| 分支 | 中位相对偏差 | 最大相对偏差 | 最大偏差位置 |
|---|---:|---:|---:|
| Ludlow16（两种 PBH 质量函数相同） | 3.18% | 13.58% | \(z=6\) |
| Prada12（两种 PBH 质量函数相同） | 1.53% | 6.97% | \(z=0\) |

剩余差异没有通过经验常数消除。较可能的来源包括：

1. 作者没有公开当时 HMFcalc 的完整版本和所有参数；
2. Ludlow16 拟合本身在高红移属于外推；
3. `MAH_corea.pdf` 的矢量反求存在有限绘图精度；
4. Prada 高峰高分支的截断实现可能仍有细微约定差异；
5. 论文内部可能使用了比图中 13 个点更密的红移网格后再绘图。

这些属于尚未完全确定的数值细节，不能伪装成已确认的物理解释。

---

## 10. 代码中的解决方法与职责边界

### 10.1 HMF

```python
press_schechter_halo_mass_function_hmfcalc(...)
```

职责：

- 公开输入、输出都使用物理 \(M_\odot\)；
- 内部转换到 `hmf` 的 \(M_\odot/h\)；
- 返回 \(dn/d\ln M\)，单位 \({\rm Mpc}^{-3}\)；
- 不负责浓度、NFW、PBH 捕获或质量吸积史。

### 10.2 Prada 精确 \(\sigma\) 分支

```python
hc.linear_sigma_hmfcalc_wmap5(...)
hc.concentration_prada12_hmf_sigma(...)
```

职责：

- 用 WMAP5 功率谱计算 \(\sigma(M,z)\)；
- 代入 Prada12 Eqs. (12)--(22)；
- 根据目标论文的处理限制高峰高上翘分支；
- 不删除或改写 `haloconcentration.py` 中的 Eq. (23) 教学实现。

### 10.3 Fig. 8 质量轨迹

下面的函数以及三组 `FIGURE8_MAH_*` 数据现位于
`notebooks/2bodycaptot.ipynb`：

```python
_figure8_mass_history(...)
```

职责：

- Ludlow16：插值论文中间矢量轨迹恢复的 \(\alpha\)、\(\beta\)；
- Prada12：由修正后的 \(C(M_0,0)\) 和 Appendix C 构造质量轨迹；
- 只服务于 Fig. 8 论文复现。

### 10.4 两套总率接口

```python
comoving_capture_rate_eq18_gpc3_per_year(...)
```

用于严格的 Eq. (18) 当前质量积分。

```python
comoving_capture_rate_with_history_gpc3_per_year(
    ...,
    mass_history_function=_figure8_mass_history,
)
```

用于论文 Fig. 8 的数值复现。总率函数保留在 `twobodycapture.py`，具体的 Fig. 8
吸积史由 Notebook 通过 `mass_history_function` 传入。函数注释中明确写出了两套
质量网格及下标配对，避免以后把它误认为一般物理公式。

`capture_rate_per_halo_a14_per_year` 新增了可选的
`concentration_override`，只是为了让 Prada 精确 \(\sigma\) 分支复用同一套 NFW、
速度分布和 Appendix A14 计算。未传入时，原有调用方式完全不变。

---

## 11. notebook 与输出文件

主复现 notebook：

```text
notebooks/2bodycaptot.ipynb
```

它完成：

1. 构造单色和对数正态 PBH 质量分布；
2. 调用主体脚本计算四条 Fig. 8 曲线；
3. 加载从原始矢量图恢复的参考点；
4. 输出端点误差和最大相对误差；
5. 保存论文风格图和诊断叠图。

输出文件：

```text
fig8_reproduction.csv
fig8_reproduction.png
fig8_reproduction_comparison.png
```

其中 `fig8_reproduction_comparison.png` 的实线是本次计算，空心圆是论文原始矢量
数据。它是检查剩余偏差最直接的文件。

整本执行命令为：

```powershell
python -m jupyter nbconvert --to notebook --execute --inplace `
    "notebooks\2bodycaptot.ipynb" `
    --ExecutePreprocessor.timeout=900
```

---

## 12. 今后应该怎样选择接口

### 情形 A：研究真实的总并合率或改变物理参数

使用：

```python
comoving_capture_rate_eq18_gpc3_per_year
```

原因：它在同一个当前晕质量上评价每晕率和 HMF，符合 Eq. (18) 的直接含义。

### 情形 B：要求尽可能贴近论文 Fig. 8

使用：

```python
comoving_capture_rate_with_history_gpc3_per_year
```

原因：它复现了从原始矢量曲线反推的轨迹和网格约定。

### 情形 C：只研究 Prada12 浓度公式

- 若学习论文闭式近似，使用 `haloconcentration.py` 中的 Eq. (23) 路径；
- 若比较目标论文 Fig. 8，使用
  `hc.concentration_prada12_hmf_sigma` 的功率谱 \(\sigma\) 路径。

---

## 13. 如果以后再次出现数量级差异，按什么顺序排查

建议固定按以下顺序，而不是直接调归一化：

1. **积分变量**：当前质量 \(M\) 还是今天的标签 \(M_0\)？
2. **质量网格身份**：每晕率和 HMF 是否真的在同一个质量点求值？
3. **HMF 单位**：\(M_\odot\) 与 \(M_\odot/h\)，\({\rm Mpc}^{-3}\) 与
   \(h^3{\rm Mpc}^{-3}\) 是否正确转换？
4. **积分测度**：使用的是 \(dn/dM\)、\(dn/d\ln M\) 还是
   \(dn/d\log_{10}M\)？
5. **质量吸积史**：固定的是同一条 \(M_0\) 轨迹，还是每个红移重新排序？
6. **浓度模型**：Prada 使用功率谱 \(\sigma\) 还是 Eq. (23) 近似？
7. **每晕率**：NFW、速度积分、PBH 双质量积分和单位是否正确？
8. **绘图数据**：画的是刚运行生成的 CSV，还是旧文件？

本次排查说明，只有完成这种分层对照，才能区分：

- 真正的物理公式错误；
- 单位或 Jacobian 错误；
- 模型近似选择不同；
- 作者未公开的数值实现约定；
- 单纯的旧文件或绘图缓存问题。

---

## 14. 最终结论

1. HMF 的 \(10^3\!-\!10^{15}M_\odot\) 上下界设置正确。
2. Bird Fig. 2 交叉检验排除了 HMF 和单晕率链中的数量级单位错误。
3. 旧 Fig. 8 的高红移差异主要来自积分对象不同：严格当前质量积分不能直接复现
   作者实际画出的轨迹/网格结果。
4. Prada12 还存在 Eq. (23) 近似使用场景错误，现已改为 WMAP5 功率谱
   \(\sigma(M,z)\) 加 Eqs. (12)--(22)。
5. Fig. 8 专用复现接口把作者式数值约定与严格 Eq. (18) 分开，避免把潜在的网格
   实现问题写成物理定律。
6. 修正后，Ludlow16 和 Prada12 四条曲线的中位偏差降到约 1.5%--3.2%，最大
   偏差为 13.6%，不再存在一个数量级的失配。

这个结果已经足以说明旧 Fig. 8 为什么失败、哪些模块是正确的、哪些近似需要
更换，以及今后应如何在“物理计算”和“论文图像复现”之间选择正确的接口。
