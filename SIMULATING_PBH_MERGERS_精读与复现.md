# 《Simulating Binary Primordial Black Hole Mergers in Dark Matter Halos》精读与复现

论文：M. Aljaf 与 I. Cholis，arXiv:2408.06515v2。论文研究暗物质晕中 PBH 双星并合，包含两条晚期通道：两体 GW 捕获，以及早期形成的硬双星在晕内经历 binary-single 相互作用。除非另述，论文取 $f_{\rm PBH}=1$，晕中 PBH 双星和单 PBH 各占一半质量；晕质量范围为 $10^3-10^{15}M_\odot$。

## 一、晕模型（Eqs. 1–10）

### 1. NFW 密度剖面

$$
\rho_{\rm NFW}(r)=\frac{\rho_s}{(r/R_s)(1+r/R_s)^2}.
$$

令 $x=r/R_s$，则 $\rho=\rho_s/[x(1+x)^2]$。因此

$$
\rho\propto r^{-1}\ (r\ll R_s),\qquad \rho\propto r^{-3}\ (r\gg R_s).
$$

晕中心密度高、外围密度快速下降。捕获需要两个 PBH，局部反应率含 $n_1n_2\propto\rho^2$，故中心区域会被强烈加权。

### 2. 浓度参数与质量归一化

浓度定义为

$$
C=\frac{R_{\rm vir}}{R_s},\qquad R_s=\frac{R_{\rm vir}}C,
$$

并定义

$$
g(C)=\ln(1+C)-\frac C{1+C}.
$$

NFW 质量积分为

$$
M(<R)=4\pi\rho_sR_s^3\left[\ln(1+R/R_s)-\frac{R/R_s}{1+R/R_s}\right].
$$

在 $R=R_{\rm vir}$ 处：

$$
M_{\rm vir}=4\pi R_s^3\rho_sg(C),\qquad
\rho_s=\frac{M_{\rm vir}}{4\pi R_s^3g(C)}.
$$

virial 半径按平均密度 $200\rho_{\rm crit}$ 定义：

$$
M_{\rm vir}=\frac{4\pi}{3}R_{\rm vir}^3(200\rho_{\rm crit}).
$$

联立 $R_{\rm vir}=CR_s$ 得

$$
\rho_s=\rho_{\rm crit}\delta_c\qquad
\delta_c=\frac{200C^3}{3g(C)}.
$$

所以 $\delta_c$ 不是独立拍脑袋的参数，而由 NFW 质量归一化和 virial 定义共同确定。

### 3. $C(M,z)$ 的物理意义

论文比较 Ludlow16 与 Prada12 两种浓度模型。一般而言，小质量晕更早形成、浓度更高；大质量晕形成较晚、浓度较低。捕获率中出现

$$
R_{\rm halo}\propto
\frac{M_{\rm vir}^2}{R_s^3g(C)^2}f(C),
\qquad f(C)=1-(1+C)^{-3}.
$$

固定总质量时，增大 $C$ 会减小 $R_s$、提高中心密度和 $\int r^2\rho^2dr$，通常增强捕获率。两种模型的差异本质上是晕集中程度的差异。

### 4. 晕质量函数

宇宙中晕质量不是单一值，需要质量函数 $dn/dM$。Press–Schechter 结构为

$$
\frac{dn}{dM}=\frac{\rho_{m,0}}M\frac{d\ln\nu}{d\ln M}f(\nu),
\qquad f(\nu)=\sqrt{\frac2\pi}\nu e^{-\nu^2/2},
$$

其中

$$
\nu(M,z)=\frac{\delta_{\rm sc}(z)}{\sigma(M,z)}.
$$

$\nu$ 大表示罕见的高密度峰，$\nu$ 小表示更常见的晕。宇宙总率为

$$
R(z)=\int dM\,R_{\rm halo}(M,z)\frac{dn}{dM}.
$$

必须区分“每晕率”和“单位共动体积率”：后者还要乘晕数密度并对质量积分。

### 5. PBH 速度分布

晕内 PBH 相对速度用截断 Maxwell 分布：

$$
P(v)=F_0v^2\left[e^{-v^2/v_{\rm disp}^2}-e^{-v_{\rm vir}^2/v_{\rm disp}^2}\right],
\qquad 0<v<v_{\rm vir},
$$

其中 $F_0$ 由

$$
4\pi\int_0^{v_{\rm vir}}v^2\left[e^{-v^2/v_{\rm disp}^2}-e^{-v_{\rm vir}^2/v_{\rm disp}^2}\right]dv=F_0^{-1}
$$

确定。论文取 $x_{\rm max}=r_{\rm max}/R_s=2.1626$：

$$
v_{\rm disp}=\sqrt{\frac{GM(<r_{\rm max})}{r_{\rm max}}},
\qquad v_{\rm vir}=\sqrt{\frac{GM_{\rm vir}}{R_{\rm vir}}},
$$
$$
M(<r_{\rm max})=4\pi R_s^3\rho_s\left[\ln(1+x_{\rm max})-\frac{x_{\rm max}}{1+x_{\rm max}}\right].
$$

由于 GW 捕获截面满足 $\sigma(v)\propto v^{-18/7}$，低速 PBH 对捕获率贡献特别大。因此存在竞争：高密度增强捕获，高速度抑制捕获；大质量晕不一定具有最大的每晕率。

## 二、两体 GW 捕获（Eqs. 11–18，Appendix A）

两质量 PBH 的捕获截面为

$$
\sigma(v)=2\pi\left(\frac{85\pi}{6\sqrt2}\right)^{2/7}
\frac{G^2(m_1+m_2)^{10/7}(m_1m_2)^{2/7}}{c^{10/7}v^{18/7}}.
$$

局部率为 $d\Gamma/d^3x=n_1n_2\langle\sigma v\rangle$。同质量时 unordered pair 因子为 $1/2$：

$$
\frac12n^2=\frac12\left[\frac{f_{\rm PBH}f_m\rho_{\rm NFW}(r)}m\right]^2.
$$

对晕体积积分：

$$
\Gamma_{\rm cap}=4\pi\int_0^{R_{\rm vir}}dr\,r^2\frac12n^2\langle\sigma v\rangle.
$$

径向部分满足

$$
\int_0^{R_{\rm vir}}r^2\rho_{\rm NFW}^2dr
=\frac{R_s^3\rho_s^2}{3}f(C),\qquad f(C)=1-(1+C)^{-3}.
$$

定义

$$
D=\int_0^{v_{\rm vir}}P(v)\left(\frac{2v}{c}\right)^{3/7}dv,
$$

则单色质量的每晕率为

$$
R_{\rm halo}=\frac{2\pi}{3}\left(\frac{85\pi}{6\sqrt2}\right)^{2/7}
\frac{G^2M_{\rm vir}^2}{R_s^3c^{10/7}g(C)^2}
f_{\rm PBH}^2f_m^2f(C)D.
$$

单色情形中 $m$ 抵消；扩展质量必须在 $5-150M_\odot$ 上积分质量核

$$
\frac{(m_1+m_2)^{10/7}(m_1m_2)^{2/7}}{m_1m_2}\phi(m_1)\phi(m_2).
$$

## 三、binary-single 动力学（Eqs. 19–32）

只演化硬双星 $a\le a_h$。

$$
\dot a=-\frac{GH\rho_{\rm env}}{v_{\rm disp}}a-
\frac{64G^3(m_1+m_2)m_1m_2}{5c^5a^3}F(e),
$$
$$
F(e)=\frac{1+73e^2/24+37e^4/96}{(1-e^2)^{7/2}},
$$
$$
\dot e=\frac{GHK\rho_{\rm env}}{v_{\rm disp}}a^{-1}-
\frac{304G^3(m_1+m_2)m_1m_2}{15c^5a^4}D(e),
\quad D(e)=\frac{e+121e^3/304}{(1-e^2)^{5/2}}.
$$

$$
\rho_{\rm env}=\rho_{\rm NFW},\quad
H=14.55(1+0.287a/a_h)^{-0.95},\quad
a_h=\frac{Gm_1}{4v_{\rm disp}^2},\quad
v_{\rm disp}=\sqrt{\frac{2GM(<r,t)}{r(t)}}.
$$

第一项是环境三体硬化/偏心率激发，第二项是 Peters GW 收缩/圆化。论文忽略 PBH 从晕中弹出，小晕率可能因此偏高。壳边界按对数划分，局部 Euler 步长为 2 Myr，全局环境更新为 200 Myr；$a\to0$ 判定并合。壳级重标定率为

$$
R_{\rm halo}=\sum_{i,t}\frac{N_{{\rm BBH},i,t}}{N_{\rm sample}}
\frac{N_{{\rm merger},i,t}}{t_{\rm look}}.
$$

内壳未重标定事件多，但真实双星数主要在外壳；重标定后外壳常主导。大质量晕因速度过高而抑制硬化。

## 四、早期双星初始分布（Appendix D）

$$
\bar x=\left(\frac{3m}{4\pi f_{\rm PBH}\rho_{\rm eq}}\right)^{1/3},\qquad j=\sqrt{1-e^2},
$$
$$
P(j)=\frac{y^2}{j(1+y^2)^{3/2}},\qquad
y=\frac{j}{0.5\sqrt{1+\sigma_{\rm eq}^2/f_{\rm PBH}^2}(x/\bar x)^3},
$$
$$
P(a,j)=\frac{3a^{-1/4}f_{\rm PBH}^{3/4}}{4\alpha\bar x}P(j)e^{-[x(a)/\bar x]^3},
\quad x(a)=\left(\frac{3am}{4\pi\alpha\rho_{\rm eq}}\right)^{1/4},\quad \alpha=0.1.
$$

高偏心率、小半长轴的双星最容易早并合；论文对 $m=30M_\odot,f=1$ 给出的 pristine 平均半长轴约为 $7287$ AU。



## 五、第一部分晕模型自检题

1. **为什么 $\rho_s$ 不是完全自由的？** 因为 NFW 质量积分必须等于 $M_{\rm vir}$，且 $R_{\rm vir}$ 按 $200\rho_{\rm crit}$ 定义。
2. **为什么捕获率含 $\rho^2$ 而非 $\rho$？** 因为一次事件需要两个 PBH，局部率为 $n_1n_2\langle\sigma v\rangle$，而 $n\propto\rho$。
3. **为什么高浓度通常提高捕获率？** 高浓度意味着更小 $R_s$ 和更高中心密度，从而增大 $\int r^2\rho^2dr$。
4. **为什么大质量晕不一定率最高？** 大质量晕虽质量大，但速度弥散通常更高；$\sigma\propto v^{-18/7}$ 会显著抑制高速捕获。

## 六、复现边界与注意事项

目录没有作者的完整 HMFcalc 表、Ludlow16/Prada12 实现、三体 $K(r,t)$ 拟合和 Monte-Carlo 原始数据。因此目前可严格复现的是：解析捕获率结构、Appendix D 初始分布、Peters/ODE 算法及图像定性趋势；要逐点重现图 1–17，还需补齐这些输入。

- Eq. (15) 的 $m$ 抵消只对单色质量成立；扩展质量必须使用 Appendix A 的质量核。
- 论文给出的 $f_{\rm PBH}^2$ 是固定环境和双星比例下的重标度，不应当视为任意丰度的普适律。
- Eq. (32) 中 $N_{\rm BBH}/N_{\rm sample}$ 是真实壳内双星数的重标定，不是简单事件数归一化。
- RVV 的早期三体形成与本文晚期 binary-single 是不同事件历史，不能仅因都叫“三体”就直接相加。

## 七、2026-09-08：Fig. 1–2 的公式链与复现记录

今天完成的代码位于 `haloconcentration.py`。这一阶段不只是画两张曲线，
而是建立后续所有晕内并合率都会重复使用的计算链：

$$
M_0
\longrightarrow C(M_0,0)
\longrightarrow (z_{-2},\alpha,\beta)
\longrightarrow M(z)
\longrightarrow C[M(z),z].
$$

因此 Fig. 2 中图例写出的 $10^3,10^6,10^9,10^{12}M_\odot$ 是今天
$z=0$ 的质量 $M_0$，不是要求晕在全部红移上保持这个质量。计算浓度时必须
先由质量吸积史得到当时的 $M(z)$。

### 7.1 Ludlow16 平滑破幂律

Ludlow16 Appendix C 将浓度写成

$$
c(\nu,z)=c_0(z)
\left(\frac{\nu}{\nu_0(z)}\right)^{-\gamma_1(z)}
\left[1+\left(\frac{\nu}{\nu_0(z)}\right)^{1/\beta(z)}\right]
^{-\beta(z)[\gamma_2(z)-\gamma_1(z)]}.
$$

这不是两个彼此独立的幂律在某一点生硬拼接。令

$$q=\frac{\nu}{\nu_0},$$

则 $q\ll1$ 时方括号趋近于 1，浓度主要由 $q^{-\gamma_1}$ 控制；
$q\gg1$ 时方括号也贡献一个幂次，使总渐近斜率变成 $-\gamma_2$。
$\beta$ 控制两种渐近行为之间转换得多快，$\nu_0$ 控制转折出现的位置，
$c_0$ 控制整条关系的大致高度。

需要特别注意：$c_0$ 不是 $\nu=\nu_0$ 时的最终浓度。此时方括号等于 2，
所以仍存在平滑转折因子。

代码首先计算

$$
\xi=\left(\frac{M}{10^{10}h^{-1}M_\odot}\right)^{-1},
$$

再依次计算 $\sigma(M,z)$、$\nu=1.686/\sigma$、五个红移参数和最终浓度。

### 7.2 Ludlow16 高红移公式的实际问题

Appendix C6 印刷出的式子为

$$
\nu_0(z)=\frac{
4.135-0.564a^{-1}-0.210a^{-2}+0.0557a^{-3}-0.00348a^{-4}}
{D(z)},\qquad a=(1+z)^{-1}.
$$

直接计算发现，分子多项式在 $z\simeq7.8$ 附近穿过零。此后 $\nu_0<0$，
而浓度公式含有 $q$ 的非整数次幂，因此会产生无效数值。这个问题不是 Python
语法错误，而是把闭式拟合外推到高红移时出现的数学病态。PBH 论文的 Fig. 1–2
却画到了 $z=12$，但具体数值处理尚未确认；仅凭这一点不能断言作者采用了
某种延拓。此处的高红移浓度问题与第 7.4.1 节的 Fig. 9 吸积参数差异需要分开。

当前复现采用透明的临时约定：

- $z\le6$ 使用出版公式；
- $z>6$ 从 $z=6$ 的 $\nu_0$ 连续接出 $\nu_0\propto(1+z)^{-1}$；
- 原始出版公式仍由单独函数保留，未被覆盖；
- 该延拓只服务于 PBH 图像复现，不能称为 Ludlow16 新公式。

这样能够复现论文中浓度在高红移继续缓慢下降、不同质量曲线逐渐靠近的趋势，
但在连接点附近仍有轻微斜率变化，需要在后续数字化比较中量化。

### 7.3 Prada12 的计算层次

Prada12 先定义时间变量

$$
x=\left(\frac{\Omega_{\Lambda0}}{\Omega_{m0}}\right)^{1/3}a.
$$

WMAP5 参数给出 $x(z=0)=1.3931$。随后通过 Eq. 12 的积分计算增长解，
并按照论文文字定义归一化为 $D(0)=1$。如果漏掉这一步，Fig. 2 中 Prada12
的浓度会整体偏低。

质量方差为

$$
\sigma(M,z)=D(z)
\frac{16.9y^{0.41}}{1+1.102y^{0.20}+6.22y^{0.333}},
\qquad
y=\left(\frac{M}{10^{12}h^{-1}M_\odot}\right)^{-1}.
$$

再由 $c_{\min}(x)$ 和 $\sigma^{-1}_{\min}(x)$ 得到 $B_0(x),B_1(x)$，最后计算

$$
\sigma'=B_1(x)\sigma(M,z),
$$

$$
\mathcal C(\sigma')=
2.881\left[\left(\frac{\sigma'}{1.257}\right)^{1.022}+1\right]
\exp\left(\frac{0.060}{\sigma'^2}\right),
$$

$$
C_{\rm Prada}=B_0(x)\mathcal C(\sigma').
$$

原始 Prada12 是 U 形关系，高峰高端会重新上翘。PBH 论文明确说明在高红移
率计算中把这一上翘分支限制到该红移的最小浓度。代码因此同时保留：

- `cap_high_peak=False`：Prada12 原始 U 形公式；
- `cap_high_peak=True`：PBH 论文使用的高峰高端截断。

Fig. 2 使用第二种口径。

### 7.4 质量吸积史

PBH 论文 Appendix C 给出

$$
z_{-2}=\left[
\frac{200C(M_0,0)^3g(1)}{A_{\rm cosmo}\Omega_{m0}g[C(M_0,0)]}
-\frac{\Omega_{\Lambda0}}{\Omega_{m0}}
\right]^{1/3}-1,
$$

$$
\beta_{\rm MAH}=-\frac{3}{1+z_{-2}},
$$

$$
\alpha_{\rm MAH}=
\frac{\ln[g(1)/g(C)]-\beta_{\rm MAH}z_{-2}}
{\ln(1+z_{-2})},
$$

其中 $A_{\rm cosmo}=798$，最终

$$
M(z)=M_0(1+z)^{\alpha_{\rm MAH}}e^{\beta_{\rm MAH}z}.
$$

代码特意使用 `beta_mah` 和 `ludlow_transition_beta` 两个名字，防止把质量吸积史
的 $\beta_{\rm MAH}$ 与 Ludlow16 平滑转折宽度 $\beta$ 混为同一个物理量。

### 7.4.1 Fig. 9 高质量端差异：证据与待证推测（2026-09-13 记录）

本节记录 2026-09-12 对原图、论文源文件和代码的只读复核结果。确认的是
**原图高质量端参数与按公开公式及当前宇宙学参数直接计算的结果不一致**；
“作者对高质量端进行了未说明的参数延拓”是有证据支持的推测，尚不是已确认的
作者实现。不能由此断言作者在某个高红移阈值切换了公式，也不能把猜测写成
Ludlow16 或 Correa15 的理论规定。

#### 来源与原图数据恢复方法

- PBH 论文：[arXiv:2408.06515 版本记录](https://arxiv.org/abs/2408.06515)，
  [v2 正文](https://arxiv.org/html/2408.06515v2)，重点为 Section II.B、Eq. (31)、
  Fig. 9 和 Appendix C。正文明确指定 Ludlow16 的 Eq. C1 用于浓度关系，
  Appendix C 指定 $A_{\rm cosmo}=798$。
- 作者提交的 [v1 源文件包](https://arxiv.org/src/2408.06515v1) 和
  [v2 源文件包](https://arxiv.org/src/2408.06515v2) 分别包含
  `figures/MAH_corea.pdf`、`MAH_corea.pdf`；LaTeX 正文将它作为 Fig. 9 插入。
  两版原始图的页面绘图内容一致，也与本地论文内嵌图一致。所检查的两个源文件包
  未包含生成该图的计算代码或数值参数表；这不等于证明作者在其他地方没有公开代码。
- 公式交叉核对：[Ludlow et al. (2016), Appendix C](https://academic.oup.com/mnras/article/460/2/1214/2608988)，
  [Correa et al. (2015), Appendix C](https://arxiv.org/html/1501.04382#A3)。

恢复对象是 PDF 保存的矢量折线顶点，而不是对截图进行像素估读。图例覆盖在曲线
之上，但底层仍保留完整路径；根据图例的颜色、线型和绘制顺序对应各个 $M_0$。
原始图内横轴 $z=1,10$ 的刻度坐标为 $x=211.165798,408.979644$，
纵轴 $M=10^3,10^{15}M_\odot$ 的刻度坐标为 $y=96.968993,314.008306$。
这些是原始 `MAH_corea.pdf` 的绘图坐标，纵坐标向上增加，不是论文页面的屏幕坐标。
映射为

$$
\log_{10}z=\frac{x-211.165798}{408.979644-211.165798},\qquad
\log_{10}(M/M_\odot)=3+
12\frac{y-96.968993}{314.008306-96.968993}.
$$

每条曲线使用 49 个有效正红移顶点，覆盖约 $z=0.244898$ 到 $12$。
排除路径最前面横坐标为 $-1$ 的对数零点/裁剪处理顶点：它不是一个可按上述映射
恢复的物理红移样本。然后在线性形式下反求参数：

$$
\log_{10}\frac{M(z)}{M_0}
=\alpha\log_{10}(1+z)+\frac{\beta z}{\ln 10}.
$$

同时放开归一化项拟合时，也能恢复图例给出的 $\log_{10}(M_0/M_\odot)$。
原始矢量顶点上的最大拟合残差约为 $4\times10^{-8}$ dex；这是绘图坐标恢复精度，
不代表模型物理精度、模拟误差或原始数据的有效位数。分别用 $z<5$、$z>5$ 的
顶点拟合，得到相同参数，说明这两段可由同一组常数描述。

#### 原图参数与当前代码的对照

代码列使用 $\Omega_{m0}=0.308$、$\Omega_{\Lambda0}=0.692$、$h=0.678$、
$A_{\rm cosmo}=798$，以及物理质量 $M_\odot$ 输入。浓度拟合内部将其转换为
以 $h^{-1}M_\odot$ 为单位的数值，即 $hM_{\rm phys}$。独立重算印刷公式得到
相同的代码列结果。

| $M_0/M_\odot$ | 原图反求 $\alpha$ | 原图反求 $\beta$ | 代码 $\alpha$ | 代码 $\beta$ |
|---:|---:|---:|---:|---:|
| $10^{13}$ | 0.270200 | -0.904000 | 0.249546 | -0.869824 |
| $10^{14}$ | 0.331200 | -1.014000 | 0.267306 | -1.087204 |
| $10^{15}$ | 0.392200 | -1.124000 | 0.235835 | -1.386694 |

三个原图质量点在恢复精度内满足 $\Delta\alpha=0.061$、$\Delta\beta=-0.110$，
可以写成下列**仅描述已提取三个点**的关系：

$$
\alpha=0.061\log_{10}(M_0/M_\odot)-0.5228,\qquad
\beta=-0.110\log_{10}(M_0/M_\odot)+0.526.
$$

这支持按质量另行指定或延拓参数的推测，但不证明作者代码实际使用了该表达式，
也不能确定其连接质量、质量点之间的插值或区间外适用性。在 $z=12$，原图这三条
曲线的 $\log_{10}(M/M_\odot)$ 分别为 $8.589761,9.084443,9.579125$，
相邻间隔约 $0.494682$ dex。当前未排序的代码轨迹则分别约为
$8.744861,8.631763,8.035906$；$10^{15}/10^{14}M_\odot$ 轨迹约在
$z=7.464$ 相交，$10^{14}/10^{13}M_\odot$ 轨迹约在 $z=10.795$ 相交。

#### 不依赖具体浓度拟合的一致性检查

将浓度 $C$ 当作未知数，先由原图 $\beta$ 和 PBH Appendix C3 得到
$z_{-2}=-3/\beta-1$，再分别用 C2、C1 反求同一个 $C$：

$$
g(C)=g(1)\exp[-\alpha\ln(1+z_{-2})-\beta z_{-2}],
$$

$$
\frac{200C^3g(1)}{A_{\rm cosmo}g(C)}
=\Omega_{m0}(1+z_{-2})^3+\Omega_{\Lambda0}.
$$

在上述当前宇宙学参数下，结果如下：

| $M_0/M_\odot$ | 原图 $\beta$ 导出的 $z_{-2}$ | C2+C3 要求的 $C$ | C1+C3 要求的 $C$ |
|---:|---:|---:|---:|
| $10^{13}$ | 2.318584 | 6.395405 | 6.587899 |
| $10^{14}$ | 1.958580 | 5.176050 | 5.740696 |
| $10^{15}$ | 1.669039 | 4.309532 | 5.080786 |

因此这些原图参数不能在该设定下同时满足 C1–C3。这个检查没有用到 $C(M)$ 的
具体拟合，单独改变质量单位换算或换一个初始浓度模型无法使同一组参数同时满足
这三个方程。若固定 $A=798$ 并保留平直宇宙条件，让 $\Omega_{m0}$ 自由取值，
三个质量点分别要求约 $0.285415,0.234238,0.195211$；也不能通过一次统一的
$\Omega_{m0}$ 替换同时消除差异。

#### 推测边界与代码处理原则

目前证据支持“高质量端存在与所述公式流程不同的参数选择”。它可能来自未说明的
经验延拓、独立参数表、实现差异或图文版本不一致；仍需作者的生成代码或说明来
确认机制。不同版本的原始图一致，只排除了本次检查的 v1/v2 之间图像变化。

高红移处的质量差异来自整条轨迹所用常数的差异，现有证据未显示某个高红移阈值
上的分段切换。第 7.2 节代码自设的 $z>6$ 浓度延拓也无法修复 Fig. 9：吸积参数
只取 $C(M_0,0)$，而不是在每个红移重新代入 $C[M(z),z]$ 求常数。

本次只补充注释和证据，保留公开公式的计算结果，并将 Fig. 9 标为尚未精确复现。
不得用逐红移排序消除交叉：那会交换不同 $M_0$ 的轨迹身份。上述反求参数也没有
被写入计算流程；若以后采用它们，应明确标为图像反演的经验输入并另行验证。

### 7.5 当前人工对照结果

当天生成：

- `fig1_reproduction(1).png`；
- `fig2_reproduction(1).png`；
- 三份相应的 CSV 曲线数据。

关键中间值为

$$
\nu_0(0)=3.41322,\qquad x_{\rm Prada}(0)=1.39311,
\qquad D_{\rm Prada}(0)=1.
$$

四个今天质量对应的 Ludlow16 浓度约为

$$
(25.00,19.78,14.38,8.71),
$$

Prada12 原始起点约为

$$
(35.80,21.68,13.28,7.92).
$$

Fig. 1 的 $10^{12}M_\odot$ 晕回溯至 $z=12$ 时质量约为
$3.25\times10^8M_\odot$，与论文图中的数量级一致。当前图已经复现主要曲线
形状和排序；精确百分比误差将在提取原图锚点后统一计算。

## 八、2026-09-13：9月9日质量吸积史阶段验收

本阶段代码位于 `halohistory.py`，并调用 `haloconcentration.py` 中已经完成的
Ludlow16 浓度和 Appendix C 质量吸积函数，没有复制第二份浓度模型。

### 8.1 已完成的计算链

对于给定的现今质量 $M_0$ 和红移数组，代码依次计算

$$
M_0\longrightarrow (z_{-2},\alpha,\beta)
\longrightarrow M(z)\longrightarrow C[M(z),z]
\longrightarrow R_{\rm vir}(z)\longrightarrow R_s(z).
$$

其中

$$
H(z)=H_0\sqrt{\Omega_{m0}(1+z)^3+\Omega_{\Lambda0}},
\qquad
\rho_{\rm crit}(z)=\frac{3H(z)^2}{8\pi G},
$$

$$
R_{\rm vir}(z)=
\left[\frac{3M(z)}{4\pi\,200\rho_{\rm crit}(z)}\right]^{1/3},
\qquad
R_s(z)=\frac{R_{\rm vir}(z)}{C[M(z),z]}.
$$

公开输入质量使用 $M_\odot$，半径输出使用 kpc。代码中的
$G=4.30091\times10^{-6}\,{\rm kpc}\,({\rm km/s})^2M_\odot^{-1}$，并把
$H_0=100h\,{\rm km\,s^{-1}\,Mpc^{-1}}$ 换成
$0.1h\,{\rm km\,s^{-1}\,kpc^{-1}}$，所以临界密度直接以
$M_\odot/{\rm kpc}^3$ 表示。

### 8.2 Fig. 9 与人工数值锚点

`make_figure9` 使用 $M_0=10^3,10^4,\ldots,10^{15}M_\odot$，保存
`fig9_reproduction.png` 和唯一一份对应的 `fig9_reproduction.csv`。CSV 保留
$z=0$，图上的对数横轴只绘制 49 个正红移点，直至 $z=12$。

以 $M_0=10^{12}M_\odot$ 为人工锚点，当前公开公式实现给出

$$
z_{-2}=3.183657,\qquad \alpha=0.223113,\qquad \beta=-0.717076,
$$

$$
C(0)=8.71031,\qquad
R_{\rm vir}(0)=210.72\ {\rm kpc},\qquad
R_s(0)=24.19\ {\rm kpc}.
$$

$R_{\rm vir}(0)$ 与论文 Table I 所列约 $211$ kpc 一致；同时
$C(0)R_s(0)=R_{\rm vir}(0)$。沿单条轨迹回溯时 $M(z)$ 总体下降，并且所有输出
保持正值。

### 8.3 阶段结论与复现边界

9月9日要求的 Appendix C 参数、$M(z)$、Fig. 9，以及
$M\rightarrow C\rightarrow R_{\rm vir}\rightarrow R_s$ 数据流已经实现并实际运行。
因此本阶段在“论文公开公式与可复用代码”层面完成。

但不能把它表述成 Fig. 9 的逐点精确复现：第 7.4.1 节已经证明，原图高质量端
反求出的 $\alpha,\beta$ 与当前宇宙学下的印刷公式 C1--C3 不能同时一致。代码保留
公开公式，不采用逐红移排序或未经来源确认的高质量参数去强行贴图。该差异作为
明确的复现误差边界保留，等待作者代码或进一步来源证据。
