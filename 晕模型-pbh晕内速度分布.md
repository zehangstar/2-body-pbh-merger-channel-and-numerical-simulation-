# 晕模型：PBH 晕内速度分布与截断 Maxwell 模型

本文配合阅读 Aljaf & Cholis 的 *Simulating Binary Primordial Black Hole Mergers in Dark Matter Halos*，重点解释论文 Section II.D、Section III 和 Section IV 中的速度模型。

需要先记住本文最重要的结论：

1. Maxwell 速度分布来自“速度的三个笛卡尔分量近似为相互独立的高斯变量”；
2. 速度分布中的 $v^2$ 来自三维速度空间的球壳体积元，不是一项新的动力学作用；
3. 截断来自有限质量晕只能束缚低于某个截止速度的 PBH；
4. 论文不是简单地把 Maxwell 分布在截止速度处砍断，而是采用在截止点连续降为零的 lowered Maxwell 形式；
5. Section III 的两体直接俘获计算采用一个全晕共享的速度分布，因此没有显式的半径 $r$；
6. 这不是“真实晕的速度分布与半径无关”，而是为了把解析积分做出来而采用的单区近似；
7. Section IV 的 binary-single 演化反而明确使用 $v_{\rm disp}^{\rm env}(r,t)$，同时依赖位置和时间。

---

## 1. 论文究竟需要什么速度？

考虑两颗 PBH，其速度矢量分别为

$$
\boldsymbol v_1,
\qquad
\boldsymbol v_2.
$$

两体引力波俘获截面依赖的不是它们各自相对于晕中心的速度，而是二者的相对速度

$$
\boxed{
v_{\rm rel}
=\left|\boldsymbol v_1-\boldsymbol v_2\right|
\equiv v_{\rm pbh}
}.
$$

必须区分以下三个概念：

| 符号 | 含义 | 是否为随机变量 |
|---|---|---|
| $\boldsymbol v_1,\boldsymbol v_2$ | 两颗 PBH 相对于晕参考系的速度矢量 | 是 |
| $v_{\rm pbh}=|\boldsymbol v_1-\boldsymbol v_2|$ | 两颗 PBH 相遇时的相对速率 | 是 |
| $v_{\rm disp}$ | 控制速度分布宽度的速度尺度 | 否，是给定晕模型后确定的参数 |

论文 Section II.D 给出的 $p(v_{\rm pbh})$，要描述的是大量 PBH 相遇事件的相对速率分布，而不是说每一颗 PBH 都以 $v_{\rm disp}$ 运动。

---

## 2. 普通 Maxwell 速率分布从哪里来？

### 2.1 从各向同性高斯速度场开始

设一个已经近似弛豫的无碰撞粒子系统没有整体漂移，并且三个速度分量近似独立：

$$
v_x,;v_y,;v_z.
$$

若每个分量都服从均值为零的高斯分布，那么三维速度矢量的联合概率密度具有形式

$$
f_{3}(\boldsymbol v)
\propto
\exp\left(-\frac{v_x^2+v_y^2+v_z^2}{v_d^2}\right)
=
\exp\left(-\frac{v^2}{v_d^2}\right),
$$

其中

$$
v=|\boldsymbol v|.
$$

参数 $v_d$ 控制分布宽度。不同文献有时把指数写成

$$
\exp\left(-\frac{v^2}{2\sigma_{1\rm D}^2}\right),
$$

这时

$$
v_d=\sqrt{2}\,\sigma_{1\rm D}.
$$

所以看到“velocity dispersion”时，不能只看名称判断它是一维标准差还是 Maxwell 指数中的宽度参数，必须看指数里有没有因子 $2$。

### 2.2 为什么速率分布里会多出 $v^2$？

$f_3(\boldsymbol v)$ 是三维速度空间中的概率密度。我们真正关心的是速率落在 $[v,v+dv]$ 内的概率。

速度空间中半径为 $v$、厚度为 $dv$ 的球壳体积是

$$
d^3v=4\pi v^2,dv.
$$

因此速率概率密度是

$$
P(v),dv
=4\pi v^2f_3(v),dv,
$$

即

$$
\boxed{
P_{\rm MB}(v)
\propto
v^2\exp\left(-\frac{v^2}{v_d^2}\right)
}.
$$

这里的 $v^2$ 只是速度空间的几何 Jacobian：速率越大，对应的速度方向球壳面积越大。

它并不表示高速粒子额外受到了一项正比于 $v^2$ 的力。

### 2.3 为什么相对速度仍然近似是 Maxwell 形式？

两颗 PBH 的相对速度矢量为

$$
\boldsymbol v_{\rm rel}
=\boldsymbol v_1-\boldsymbol v_2.
$$

如果 $\boldsymbol v_1$ 和 $\boldsymbol v_2$ 的每个分量都是相互独立的高斯变量，那么两个高斯变量之差仍然是高斯变量。因此 $\boldsymbol v_{\rm rel}$ 的三个分量仍为高斯分布，其模长仍具有 Maxwell 型速率分布。

若每颗 PBH 的单个速度分量方差都是 $s^2$，并且二者统计独立，则

$$
{\rm Var}(v_{1x}-v_{2x})
=s^2+s^2
=2s^2.
$$

也就是说，相对速度分布比单粒子速度分布更宽。论文直接用 $v_{\rm disp}$ 参数化相对速度分布；复现时应忠实使用论文的这一参数定义，不要再次手动乘一个 $\sqrt2$，否则可能重复计算相对速度展宽。

---

## 3. 为什么必须截断 Maxwell 分布？

### 3.1 无限 Maxwell 尾部不适合有限引力势阱

普通 Maxwell 分布定义在

$$
0\leq v<\infty.
$$

但暗物质晕是一个有限深度的引力势阱。速度足够大的粒子具有正的总能量，会逃出晕，不能继续作为稳定的晕成员。

因此束缚 PBH 的速率分布需要一个截止速度 $v_{\rm cut}$：

$$
P(v)=0,
\qquad
v>v_{\rm cut}.
$$

目标论文把这个截止速度记为 $v_{\rm vir}$。为了避免它与通常的圆周 virial velocity 混淆，理解时可以先把它读成

$$
v_{\rm cut}\equiv v_{\rm vir}^{\rm(paper)}.
$$

### 3.2 硬截断与 lowered Maxwell 的区别

最简单的硬截断是

$$
P_{\rm hard}(v)
\propto
v^2e^{-v^2/v_{\rm disp}^2}
\Theta(v_{\rm cut}-v).
$$

它在截止点内侧仍为有限正值，然后突然跳到零。

目标论文使用的是更平滑的 lowered Maxwell 形式。先在三维速度概率密度层面定义

$$
\Phi(v)=
\begin{cases}
\displaystyle
\exp\left(-\frac{v^2}{v_{\rm disp}^2}\right)
-
\exp\left(-\frac{v_{\rm cut}^2}{v_{\rm disp}^2}\right),
&0\leq v\leq v_{\rm cut},\\[8pt]
0,&v>v_{\rm cut}.
\end{cases}
$$

因为

$$
\Phi(v_{\rm cut})=0,
$$

所以概率密度会在截止点连续降到零，而不是突然被砍断。

论文 Eq. (7) 的核心结构正是

$$
\boxed{
v_{\rm pbh}^2
\left[
e^{-v_{\rm pbh}^2/v_{\rm disp}^2}
-e^{-v_{\rm vir}^2/v_{\rm disp}^2}
\right]
},
\qquad
0\leq v_{\rm pbh}\leq v_{\rm vir}.
$$

其中：

- $v_{\rm pbh}^2$ 是三维速度空间球壳因子；
- 第一项是 Maxwell 型指数衰减；
- 第二项把整个分布向下平移，使其在 $v_{\rm vir}$ 处变为零；
- 超过 $v_{\rm vir}$ 的速度被排除，因为对应 PBH 不再视为束缚于该晕。

---

## 4. 怎样无歧义地归一化？

论文 Eq. (7)-(8) 在 $p(v)$ 是“三维速度密度”还是“一维速率密度”之间存在容易混淆的记号。代码中最安全的做法是明确区分二者。

目标论文印出的 Eq. (7) 把 $v^2$ 写进了 $p(v)$，看起来像一维速率密度；但 Eq. (8) 的归一化常数又保留了三维球壳因子 $4\pi v^2$。它所引用的 Bird et al. 原式则把 $P(v)$ 写成不含 $v^2$ 的三维各向同性速度密度，再用

$$
4\pi\int P(v)v^2,dv=1
$$

归一化。因此，不应只机械复制字母 $F_0$；应先明确程序返回的是 $f_3(v)$ 还是速率 PDF $P(v)$，再决定 $4\pi v^2$ 放在哪里。

### 4.1 三维各向同性速度密度

定义

$$
f_3(\boldsymbol v)
=A\Phi(v).
$$

三维归一化要求

$$
\int f_3(\boldsymbol v),d^3v
=4\pi\int_0^{v_{\rm cut}}v^2f_3(v),dv
=1.
$$

因此

$$
\boxed{
A^{-1}
=4\pi\int_0^{v_{\rm cut}}
v^2\Phi(v),dv
}.
$$

### 4.2 一维速率概率密度

与之对应的、可以直接对 $dv$ 积分的速率概率密度是

$$
\boxed{
P(v)=4\pi v^2f_3(v)
=4\pi A v^2\Phi(v)
}.
$$

它满足

$$
\boxed{
\int_0^{v_{\rm cut}}P(v),dv=1
}.
$$

后续计算

$$
\langle Q(v)\rangle
$$

时，只能选择以下两种完全等价的写法之一：

$$
\langle Q\rangle
=\int d^3v,f_3(\boldsymbol v)Q(v),
$$

或

$$
\langle Q\rangle
=\int_0^{v_{\rm cut}}dv,P(v)Q(v).
$$

不能在 $P(v)$ 中已经包含 $4\pi v^2$ 后，又在积分外再次乘 $4\pi v^2$；也不能两个地方都不包含它。

### 4.3 归一化积分的闭式结果

令

$$
q=\frac{v_{\rm cut}}{v_{\rm disp}}.
$$

则

$$
\int_0^{v_{\rm cut}}
v^2\Phi(v),dv
=v_{\rm disp}^3
\left[
\frac{\sqrt\pi}{4}{\rm erf}(q)
-e^{-q^2}
\left(
\frac q2+\frac{q^3}{3}
\right)
\right].
$$

所以

$$
A^{-1}
=4\pi v_{\rm disp}^3
\left[
\frac{\sqrt\pi}{4}{\rm erf}(q)
-e^{-q^2}
\left(
\frac q2+\frac{q^3}{3}
\right)
\right].
$$

这可以用于以后检查数值积分，但本文不添加任何自动自检代码。

---

## 5. $v_{\rm disp}$ 如何由 NFW 晕确定？

### 5.1 NFW 晕的内部质量

NFW 密度为

$$
\rho_{\rm NFW}(r)
=\frac{\rho_s}
{(r/R_s)(1+r/R_s)^2}.
$$

定义

$$
x=\frac r{R_s},
\qquad
g(x)=\ln(1+x)-\frac{x}{1+x}.
$$

则半径 $r=xR_s$ 内的质量为

$$
\boxed{
M(<r)=4\pi\rho_sR_s^3g(x)
}.
$$

整个晕满足

$$
M_{\rm vir}
=4\pi\rho_sR_s^3g(C),
\qquad
C=\frac{R_{\rm vir}}{R_s}.
$$

因此也可以写成

$$
M(<r)
=M_{\rm vir}\frac{g(x)}{g(C)}.
$$

### 5.2 为什么选 $r_{\max}=2.1626R_s$？

圆周速度定义为

$$
V_c^2(r)=\frac{GM(<r)}r.
$$

对 NFW 晕，

$$
V_c^2(x)
\propto
\frac{g(x)}x.
$$

令 $d[g(x)/x]/dx=0$，得到最大圆周速度所在位置

$$
\boxed{
x_{\max}
=\frac{r_{\max}}{R_s}
\simeq2.1626
}.
$$

目标论文把晕的特征速度弥散取为该位置的圆周速度：

$$
\boxed{
v_{\rm disp}
=\sqrt{rac{GM(<r_{\max})}{r_{\max}}}
}.
$$

代入 NFW 内部质量，得到

$$
v_{\rm disp}^2
=rac{GM_{\rm vir}}{R_{\rm vir}}
\frac{C}{x_{\max}}
\frac{g(x_{\max})}{g(C)}.
$$

这一步非常关键：论文没有在每个半径分别求一个局部速度弥散，而是选取 $r_{\max}$ 处的一个特征速度，代表整个晕的相对速度分布宽度。

### 5.3 $v_{\rm vir}$ 的命名约定

论文 Eq. (9) 写成

$$
v_{\rm disp}
=\frac{v_{\rm vir}}{\sqrt2}
\sqrt{
\frac{C}{x_{\max}}
\frac{g(x_{\max})}{g(C)}
}.
$$

与上一式比较可知，这里的 $v_{\rm vir}$ 满足

$$
\boxed{
v_{\rm vir}^{\rm(paper)}
=\sqrt{\frac{2GM_{\rm vir}}{R_{\rm vir}}}
}.
$$

而很多晕模型文献把圆周 virial velocity 定义为

$$
V_{\rm vir}^{\rm(circ)}
=\sqrt{\frac{GM_{\rm vir}}{R_{\rm vir}}}.
$$

两者相差 $\sqrt2$：

$$
v_{\rm vir}^{\rm(paper)}
=\sqrt2,V_{\rm vir}^{\rm(circ)}.
$$

因此在代码中最好采用更明确的变量名，例如 `v_cut` 和 `v_circ_vir`，不要仅凭 `v_vir` 这个名称决定是否包含 $\sqrt2$。

---

## 6. 为什么 Section III 的速度分布与半径无关？

### 6.1 直接原因：论文只给每个晕指定一个 $v_{\rm disp}$

Section III 的截断 Maxwell 分布可概念性地写成

$$
P(v_{\rm pbh}\mid M,C,z).
$$

它通过以下全局晕参数确定：

$$
M,quad C,quad R_s,quad R_{\rm vir},quad
v_{\rm disp}(r_{\max}),quad v_{\rm cut}.
$$

其中没有把正在发生相遇的位置 $r$ 作为条件，因此

$$
\boxed{
P(v_{\rm pbh}\mid r,M,C,z)
\;longrightarrow\;
P(v_{\rm pbh}\mid M,C,z)
}.
$$

这就是它在 Section III 中与径向位置无关的直接原因：作者主动使用一个代表整个晕的单区速度模型。

这不是从 NFW 密度剖面严格推导出的结论。

### 6.2 这种近似带来了什么数学便利？

两体直接俘获的局部率为

$$
\frac{d\Gamma_{\rm cap}}{d^3x}
=n_1(r)n_2(r)
\left\langle
\Sigma_{\rm cap}(v_{\rm pbh})v_{\rm pbh}
\right\rangle.
$$

这里用 $\Sigma_{\rm cap}$ 表示俘获截面，以免与宇宙学质量方差 $\sigma(M,z)$ 混淆。

如果速度分布对整个晕相同，那么

$$
\left\langle\Sigma_{\rm cap}v\right\rangle
$$

不依赖 $r$，可以移到径向积分外面：

$$
\Gamma_{\rm cap}
=\left\langle\Sigma_{\rm cap}v\right\rangle
4\pi\int_0^{R_{\rm vir}}
r^2n_1(r)n_2(r),dr.
$$

若 PBH 数密度跟随 NFW 密度，

$$
n_i(r)\propto\rho_{\rm NFW}(r),
$$

就得到

$$
\Gamma_{\rm cap}
\propto
\left\langle\Sigma_{\rm cap}v\right\rangle
\int_0^{R_{\rm vir}}r^2\rho_{\rm NFW}^2(r),dr.
$$

于是问题被拆成两个互相独立的积分：

$$
\boxed{
\text{速度平均}
\times
\text{空间密度平方积分}
}.
$$

这正是采用全局速度分布的主要计算便利。

### 6.3 真实晕中速度分布一般会依赖半径

对于稳态、球对称的无碰撞系统，局部径向速度弥散应满足 Jeans 方程

$$
\boxed{
\frac{d[\rho(r)\sigma_r^2(r)]}{dr}
+\frac{2\beta_{\rm ani}(r)}{r}
\rho(r)\sigma_r^2(r)
=-\rho(r)\frac{GM(<r)}{r^2}
},
$$

其中 $\beta_{\rm ani}$ 是速度各向异性参数。

因为 NFW 的 $\rho(r)$ 和 $M(<r)$ 都随半径变化，所以 Jeans 方程的解通常也是

$$
\sigma_r=\sigma_r(r).
$$

此外，局部逃逸速度取决于引力势：

$$
v_{\rm esc}(r)
=\sqrt{2[\Phi(\infty)-\Phi(r)]},
$$

也通常随 $r$ 改变。

因此更真实的模型应当写成

$$
P(v_{\rm pbh}\mid r,M,C,z),
$$

并计算

$$
\Gamma_{\rm cap}
=4\pi\int_0^{R_{\rm vir}}dr,r^2n_1(r)n_2(r)
\int dv,P(v\mid r)\Sigma_{\rm cap}(v)v.
$$

这时速度积分不能再简单地从径向积分中提出去。

所以最准确的表述是：

$$
\boxed{
\text{论文的 Section III 假定速度分布与 }r\text{ 无关；真实 NFW 晕并非如此。}
}
$$

---

## 7. 单区 Maxwell 近似与 Jeans 方程的关系

### 7.1 先给出精度判断

把整个晕的相对速度分布近似为同一个截断 Maxwell 分布，并令其宽度由

$$
v_{\rm disp}=V_c(r_{\max})=V_{\max}
$$

确定，可以作为解析计算和数量级估计的基准，但不能视为高精度的局部相空间模型。

它适合：

- 快速估算单个晕的两体俘获率；
- 比较不同晕质量、浓度和红移造成的主要趋势；
- 与采用相同假设的早期工作进行比较。

它的局限是：

- 忽略真实速度弥散的径向变化；
- 忽略速度各向异性；
- 把完整速度分布强制近似为 Gaussian/Maxwell 形状；
- 忽略密度最高区域和局部低速矩之间的相关性；
- 对快速吸积、未弛豫或有子结构的晕可能不准确。

因此最合适的评价是

$$
\boxed{
\text{合理的解析基准和数量级近似，但不是精确的局部动力学模型。}
}
$$

### 7.2 它不是 Jeans 方程的正式“零阶解”

如果把真实的径向函数

$$
\sigma_r(r)
$$

直接替换成一个常数

$$
\sigma_r(r)\approx\sigma_{\rm char},
$$

那么可以非正式地称其为“对径向结构的零阶粗粒化”或“单区近似”。

但是论文并没有进行如下 Taylor 展开：

$$
\sigma_r(r)
=\sigma_0+\sigma_1(r-r_0)+\cdots,
$$

也没有从 Jeans 方程证明

$$
\sigma_0=V_c(r_{\max}).
$$

所以这里的“零阶”只能表示忽略空间变化，不能理解成一个具有严格小参数和误差阶数的微扰展开。

还要注意，“Jeans 方程是一阶速度矩方程”中的“一阶”，与“近似到一阶”不是同一个概念。前者说的是对无碰撞 Boltzmann 方程乘速度并在速度空间积分，后者才是通常所说的展开精度。

### 7.3 各向同性 NFW 晕的 Jeans 解是什么样？

球对称、稳态 Jeans 方程为

$$
\frac{d[\rho(r)\sigma_r^2(r)]}{dr}
+\frac{2\beta_{\rm ani}(r)}{r}
\rho(r)\sigma_r^2(r)
=-\rho(r)\frac{GM(<r)}{r^2}.
$$

若进一步假设各向同性，

$$
\beta_{\rm ani}(r)=0,
$$

则

$$
\frac{d[\rho(r)\sigma_r^2(r)]}{dr}
=-\rho(r)\frac{GM(<r)}{r^2}.
$$

给定外边界 $R_{\rm out}$，并采用

$$
\rho(R_{\rm out})\sigma_r^2(R_{\rm out})=0
$$

这类边界条件，可以积分得到

$$
\boxed{
\sigma_r^2(r)
=\frac{1}{\rho(r)}
\int_r^{R_{\rm out}}
\rho(s)\frac{GM(<s)}{s^2}\,ds
}.
$$

NFW 晕的 $\rho(r)$ 和 $M(<r)$ 都随半径变化，因此 Jeans 解必然一般也是

$$
\sigma_r=\sigma_r(r),
$$

而不是一个全晕常数。

改变速度各向异性 $\beta_{\rm ani}(r)$、晕的外边界条件或截断方法，还会改变解出的速度弥散剖面。

### 7.4 Jeans 方程并不能唯一给出 Maxwell 分布

Jeans 方程主要约束局部二阶速度矩

$$
\sigma_r^2(r)=\left\langle v_r^2\right\rangle.
$$

但只知道方差，不能唯一确定完整的概率分布。Gaussian、平顶分布、长尾分布以及各向异性分布都可能具有相同的二阶矩。

因此，从 Jeans 方程到局部 Maxwell 模型还需要额外假设：

$$
\boxed{
\text{Jeans 给出二阶矩}
+
\text{局部各向同性、高斯闭合}
\Longrightarrow
\text{局部 Maxwell 近似}
}.
$$

若要在球对称、各向同性和平衡条件下从密度和势求完整分布，可以使用 Eddington 反演得到

$$
f(\mathcal E),
$$

然后构造局部速率分布

$$
\boxed{
P(v\mid r)
=\frac{4\pi v^2}{\rho(r)}
f\left[\Psi(r)-\frac{v^2}{2}\right]
},
$$

其自然截止速度为

$$
v_{\rm esc}(r)=\sqrt{2\Psi(r)}.
$$

这比“Jeans 方差加局部 Maxwell”更自洽，因为它同时确定分布形状和局部截止；但它仍依赖球对称、平衡和各向同性假设。

### 7.5 为什么 $V_{\max}$ 至少是合理的数量级？

virial equilibrium 给出的典型尺度是

$$
\sigma^2\sim\frac{GM}{R},
$$

而圆周速度满足

$$
V_c^2(r)=\frac{GM(<r)}{r}.
$$

所以

$$
\sigma\sim V_c
$$

在数量级上是合理的。选择 $V_{\max}$ 还有几个实际优点：

- 它是 NFW 剖面天然定义的稳定特征速度；
- 它同时包含晕质量和浓度的信息；
- 它比任意选择一个半径的局部速度更容易跨晕比较；
- 它便于和模拟中的 $V_{\max}$ 数据及早期 PBH 文献衔接。

但比例系数并不普适。例如理想奇异等温球满足

$$
V_c^2=2\sigma_{1{\rm D}}^2.
$$

此外，论文指数中的 $v_{\rm disp}$ 是相对速度分布的宽度参数，不必严格等于单粒子的一维标准差。因而

$$
v_{\rm disp}=V_{\max}
$$

应理解成特征速度的匹配约定，而不是普适动力学恒等式。

### 7.6 为什么它对 PBH 俘获率可能不够精确？

若保留局部速度分布，单晕俘获率中的关键结构是

$$
\Gamma_{\rm exact}
\propto
\int_0^{R_{\rm vir}}
4\pi r^2\rho^2(r)
\left\langle v^{-11/7}\right\rangle_r
dr.
$$

其中：

- $\rho^2(r)$ 强烈提高中心高密度区域的权重；
- $\langle v^{-11/7}\rangle_r$ 强烈提高低相对速度事件的权重。

全晕单区模型把它近似为

$$
\Gamma_{\rm 1zone}
\propto
\left\langle v^{-11/7}\right\rangle_{\rm char}
\int_0^{R_{\rm vir}}
4\pi r^2\rho^2(r)\,dr.
$$

它实际上假设速度矩与半径无关，从而忽略了

$$
\rho^2(r)
\quad\text{与}\quad
\left\langle v^{-11/7}\right\rangle_r
$$

之间的空间相关性。

可以把局部处理相对于单区处理的修正定义为

$$
\boxed{
B_v=
\frac{
\displaystyle
\int dr\,r^2\rho^2(r)
\left\langle v^{-11/7}\right\rangle_r
}{
\displaystyle
\left\langle v^{-11/7}\right\rangle_{\rm char}
\int dr\,r^2\rho^2(r)
}
}.
$$

单区近似等价于假设

$$
B_v\simeq1,
$$

但目标论文并没有证明它对所有晕质量、红移和浓度都成立。

### 7.7 从基准模型到高精度模型的层次

| 层次 | 速度处理 | 能得到什么 |
|---|---|---|
| 1. 全晕单区 Maxwell | $P(v\mid r)\rightarrow P(v\mid V_{\max})$ | 快速解析基准 |
| 2. Jeans + 局部 Maxwell | 先求 $\sigma_r(r)$，再假设每个半径局部 Maxwell | 保留主要径向速度梯度 |
| 3. Eddington 反演 | 从 $\rho(r)$ 和 $\Phi(r)$ 求 $f(\mathcal E)$ | 各向同性平衡条件下的完整局部分布 |
| 4. 模拟相空间分布 | 直接测量 $P(v_{\rm rel}\mid r,M,z)$ | 可包含各向异性、非高斯和非平衡效应 |

因此，更准确的总结是

$$
\boxed{
\begin{aligned}
&\text{论文的 }V_{\max}+\text{截断 Maxwell 模型}\\
&\text{可以视为径向粗粒化的单区基准，}\\
&\text{但不是 Jeans 方程的正式零阶或一阶解。}
\end{aligned}
}
$$

如果完成论文原始复现后再提高精度，最自然的第一步是：

$$
\text{解各向同性 NFW Jeans 方程得到 }\sigma_r(r)
\longrightarrow
\text{构造局部相对速度分布}
\longrightarrow
\text{重新进行联合的 }(r,v)\text{ 积分}.
$$

---

## 8. 截断速度分布如何进入两体俘获率？

### 8.1 引力波俘获截面对速度非常敏感

目标论文 Eq. (11) 给出

$$
\Sigma_{\rm cap}(v_{\rm pbh})
\propto
v_{\rm pbh}^{-18/7}.
$$

一次相遇的率还要乘相对速度，所以

$$
\Sigma_{\rm cap}(v)v
\propto
v^{-11/7}.
$$

速度平均为

$$
\boxed{
\left\langle\Sigma_{\rm cap}v\right\rangle
=\int_0^{v_{\rm cut}}
P(v)\Sigma_{\rm cap}(v)v,dv
}.
$$

低速相遇拥有更大的俘获截面，因此速度分布的低速部分非常重要。

### 8.2 为什么 $v\to0$ 时积分没有真的发散？

当 $v\to0$ 时，lowered Maxwell 速率密度满足

$$
P(v)\propto v^2.
$$

因此完整被积函数的低速行为是

$$
P(v)\Sigma_{\rm cap}(v)v
\propto
v^2v^{-18/7}v
=v^{3/7}.
$$

由于

$$
\int_0^\epsilon v^{3/7}dv
$$

有限，所以速度趋近零时不会产生数学发散。

低速事件很重要，但相空间体积因子 $v^2$ 抑制了速度恰好为零的事件。

### 8.3 论文中的 $D(v_{\rm pbh})$

论文把速度平均中的相关部分整理成

$$
D(v_{\rm pbh})
=\int_0^{v_{\rm vir}}
P(v_{\rm pbh},v_{\rm disp})
\left(\frac{2v_{\rm pbh}}c\right)^{3/7}
dv_{\rm pbh}.
$$

它最终进入单晕俘获率

$$
R_{\rm halo}(M,z).
$$

这里的 $P(v_{\rm pbh},v_{\rm disp})$ 沿用了论文自己的归一化记号。真正编程时，若采用本文定义的归一化速率 PDF，就应优先直接计算

$$
\left\langle\Sigma_{\rm cap}v\right\rangle
=\int_0^{v_{\rm cut}}
P(v)\Sigma_{\rm cap}(v)v,dv,
$$

然后再与论文整理后的常数因子核对。这样比同时照搬 Eq. (7)、Eq. (8) 和 Eq. (17) 的 $F_0/P$ 记号更不容易重复计入速度空间 Jacobian。

完整逻辑为

$$
(M,z,C)
\longrightarrow
(R_s,R_{\rm vir},\rho_s)
\longrightarrow
(v_{\rm disp},v_{\rm cut})
\longrightarrow
P(v)
\longrightarrow
\langle\Sigma_{\rm cap}v\rangle
\longrightarrow
R_{\rm halo}(M,z).
$$

---

## 9. Section IV 为什么又出现了位置相关的速度？

论文在 binary-single 相互作用部分不再采用完全相同的单区处理，而是明确写出

$$
\boxed{
v_{\rm disp}^{\rm env}(r,t)
=\sqrt{\frac{2GM(<r,t)}{r}}
}.
$$

这个速度同时依赖：

- 半径 $r$，因为 $M(<r,t)$ 随位置变化；
- 时间 $t$，因为晕的质量、浓度和尺度半径会演化。

它进入 binary-single 硬化方程，例如

$$
\frac{da}{dt}
\supset
-\frac{G H(r,t)\rho_{\rm env}(r,t)}
{v_{\rm disp}^{\rm env}(r,t)}a^2,
$$

以及硬双星判据中的临界半长轴

$$
a_h(r,t)
=\frac{Gm_1}
{4[v_{\rm disp}^{\rm env}(r,t)]^2}.
$$

为保留这种位置依赖，论文把较大质量的晕划分成多个球壳，并在每个球壳中分别更新：

$$
\rho_i(t),
\qquad
v_{{\rm disp},i}(t),
\qquad
a_{h,i}(t).
$$

因此整篇论文实际上是：

| 计算通道 | 速度模型 | 是否依赖半径 |
|---|---|---|
| Section III：两体直接俘获 | 由 $r_{\max}$ 定义的全晕截断 Maxwell 分布 | 否，是单区近似 |
| Section IV：binary-single 演化 | $v_{\rm disp}^{\rm env}(r,t)=\sqrt{2GM(<r,t)/r}$ | 是，采用球壳处理 |

这也解释了为什么论文引言会强调较大晕中存在速度梯度：该梯度主要在 Section IV 的分层 binary-single 数值演化中被显式保留，而不是在 Section III 的解析两体俘获积分中完整保留。

---

## 10. 为什么最小质量晕可以只用一个区域？

论文指出，对 $O(10^3)M_\odot$ 的小晕只使用一个球形区域；随着晕质量增加，才逐渐增加球壳数，质量大于约 $10^7M_\odot$ 时使用十个球壳。

作者给出的理由是：对最小质量晕，PBH 穿越整个晕的时间小于数值模拟的时间步长，因此他们不解析稳定的径向速度梯度，而把小晕看作一个混合较快的单区系统。

这里仍应区分：

$$
\boxed{
\text{模拟时间分辨率内不追踪梯度}
\neq
\text{真实物理系统严格没有梯度}
}.
$$

它是一种与时间步长和空间分箱共同有关的数值近似。

---

## 11. 这种速度近似会带来什么影响？

### 11.1 密度平方使中心区域权重很大

两体俘获率含有

$$
n_1(r)n_2(r)
\propto
\rho_{\rm NFW}^2(r).
$$

所以中心高密度区域对径向积分非常重要。

### 11.2 速度又以负幂影响俘获

同时

$$
\Sigma_{\rm cap}v
\propto v^{-11/7}.
$$

如果中心和外部的真实速度分布不同，用一个统一的 $v_{\rm disp}$ 就会把“密度最高的位置”和“实际局部速度”之间的相关性丢掉。

更一般地，真实率应包含

$$
\rho^2(r)
\left\langle v^{-11/7}\right\rangle_r,
$$

而单区模型近似为

$$
\rho^2(r)
\left\langle v^{-11/7}\right\rangle_{\rm halo}.
$$

二者不一定相同。

### 11.3 Maxwell 形式本身也只是近似

真实 CDM 晕可能具有：

- 速度各向异性；
- 径向流入；
- 子结构；
- 并合遗留的非高斯尾部；
- 随半径变化的逃逸速度；
- PBH 离散性和质量偏析。

因此 Maxwell 模型的作用是提供一个可计算、可与早期文献比较的基准，而不是宣称真实晕的六维相空间分布严格服从热平衡 Maxwell-Boltzmann 分布。

---

## 12. 复现代码时最容易犯的错误

### 12.1 把单粒子速度当成相对速度

俘获截面需要

$$
v_{\rm pbh}=|\boldsymbol v_1-\boldsymbol v_2|,
$$

不是任意一颗 PBH 的速度大小。

### 12.2 重复加入 $v^2$ 或 $4\pi$

如果函数返回的是速率 PDF $P(v)$，应满足

$$
\int P(v)dv=1,
$$

后续积分不再乘 $4\pi v^2$。

如果函数返回的是三维各向同性密度 $f_3(v)$，则应使用

$$
4\pi\int v^2 f_3(v)dv=1.
$$

两套记号只能选一套。

### 12.3 把 $v_{\rm vir}$ 与圆周 virial velocity 混为一谈

目标论文 Eq. (9) 的约定包含一个 $1/\sqrt2$，表明其截止量满足

$$
v_{\rm cut}=\sqrt{2GM_{\rm vir}/R_{\rm vir}}.
$$

若代码已经使用这个截止速度，就不能再额外乘一次 $\sqrt2$。

### 12.4 把 $v_{\rm disp}$ 当成处处相同的物理定律

Section III 中它是全晕近似；Section IV 中则应换成局部的 $v_{\rm disp}^{\rm env}(r,t)$。

### 12.5 混淆三个不同的 $\sigma$

本项目中至少会遇到：

| 记号 | 含义 |
|---|---|
| $\sigma(M,z)$ | 线性密度涨落的 rms，用于峰高和浓度模型 |
| $\Sigma_{\rm cap}(v)$ | 两颗 PBH 的引力波俘获截面 |
| $\sigma_r(r)$ 或 $v_{\rm disp}$ | 晕内速度弥散 |

它们在物理上完全不同。代码中最好分别命名为 `sigma_density`、`capture_cross_section` 和 `velocity_dispersion`。

---

## 13. 从晕参数到两体俘获率的计算链

给定晕在红移 $z$ 时的质量和浓度：

1. 计算
   $$R_{\rm vir},\quad R_s=R_{\rm vir}/C;$$
2. 计算
   $$\rho_s,\quad M(<r);$$
3. 取
   $$r_{\max}=2.1626R_s;$$
4. 计算全晕特征速度
   $$v_{\rm disp}=\sqrt{GM(<r_{\max})/r_{\max}};$$
5. 按论文约定计算截止速度
   $$v_{\rm cut}=v_{\rm vir}^{\rm(paper)};$$
6. 构造并归一化 lowered Maxwell 速率分布 $P(v)$；
7. 计算
   $$\langle\Sigma_{\rm cap}v\rangle;$$
8. 计算空间积分
   $$4\pi\int r^2n_1(r)n_2(r)dr;$$
9. 两部分相乘，得到
   $$R_{\rm halo}(M,z).$$

流程图为

$$
(M,z,C)
\longrightarrow
(R_{\rm vir},R_s,\rho_s)
\longrightarrow
(v_{\rm disp},v_{\rm cut})
\longrightarrow
P(v)
\longrightarrow
\langle\Sigma_{\rm cap}v\rangle
\longrightarrow
R_{\rm halo}.
$$

---

## 14. 最终物理图像

截断 Maxwell 模型可以理解为以下近似链：

$$
\text{近似弛豫、各向同性的随机速度}
\longrightarrow
\text{高斯速度分量}
\longrightarrow
\text{Maxwell 速率分布}
\longrightarrow
\text{有限晕势阱的速度截断}.
$$

在目标论文的两体俘获通道中，又进一步采用

$$
\text{局部速度分布}
\longrightarrow
\text{由 }r_{\max}\text{ 定义的全晕代表速度分布},
$$

从而把速度积分与 NFW 径向积分分离。

所以，对“为什么速度分布与径向位置无关”的最终回答是：

> 它并不是物理上必然与半径无关。论文在两体直接俘获的解析计算中，用
> $r_{\max}=2.1626R_s$ 处的特征速度代表整个晕，因此人为消除了显式的
> $r$ 依赖；在 binary-single 数值演化中，作者重新引入了随半径和时间变化的
> $v_{\rm disp}^{\rm env}(r,t)$。

---

## 参考文献

1. Aljaf & Cholis, *Simulating Binary Primordial Black Hole Mergers in Dark Matter Halos*, [arXiv:2408.06515](https://arxiv.org/abs/2408.06515)，尤其是 Eqs. (7)-(10)、(12)-(17)、(19)、(23)-(25)。
2. Bird et al., *Did LIGO Detect Dark Matter?*, [Phys. Rev. Lett. 116, 201301 (2016)](https://arxiv.org/abs/1603.00464)，尤其是其截断 Maxwell 速度模型与两体俘获率 Eqs. (5)-(9)。
3. Cholis et al., *Orbital eccentricities in primordial black hole binaries*, [Phys. Rev. D 94, 084013 (2016)](https://arxiv.org/abs/1606.07437).
4. Mao et al., *Halo-to-Halo Similarity and Scatter in the Velocity Distribution of Dark Matter*, [Astrophys. J. 764, 35 (2013)](https://arxiv.org/abs/1210.2721)，讨论模拟晕速度分布相对于简单 Maxwell 模型的偏离。
5. Łokas & Mamon, *Properties of spherical galaxies and clusters with an NFW density profile*, [MNRAS 321, 155 (2001)](https://arxiv.org/abs/astro-ph/0002395)，给出不同各向异性假设下 NFW 晕的径向速度弥散解。
6. Wojtak et al., *Radial velocity moments of dark matter haloes*, [MNRAS 361, L1 (2005)](https://arxiv.org/abs/astro-ph/0503391)，比较模拟晕的径向速度矩、非高斯性与 Jeans 解。
