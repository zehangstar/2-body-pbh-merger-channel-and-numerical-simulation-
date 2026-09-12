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
却画到了 $z=12$，但没有说明作者如何延拓该式。

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

### 7.5 当前人工对照结果

今天生成：

- `fig1_reproduction.png`；
- `fig2_reproduction.png`；
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
