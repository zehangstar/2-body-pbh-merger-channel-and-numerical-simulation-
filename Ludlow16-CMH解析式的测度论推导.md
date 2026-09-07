# 从高斯随机场到 Ludlow16 的 collapsed mass history

这份教程用测度论和随机过程语言解释 Ludlow et al. (2016) 中的解析式

$$
\frac{M_{\rm coll}(z;f)}{M_0}
=
\operatorname{erfc}\!\left[
\frac{\delta_{\rm sc}(z)-\delta_{\rm sc}(z_0)}
{\sqrt{2\left[\sigma^2(fM_0)-\sigma^2(M_0)\right]}}
\right].
$$

本文沿用《随机过程-高斯随机场-Euclidean场论教程》的基本观念：

> 概率密度不是最基本的对象，概率测度才是；随机过程是一条路径值随机变量。

我们要建立的完整链条是

$$
\boxed{
(\Omega,\mathcal F,P)
\longrightarrow
\text{Gaussian random field}
\longrightarrow
\delta_R(\mathbf q)
\longrightarrow
\delta(S)
\longrightarrow
\text{Brownian motion}
\longrightarrow
\text{first crossing}
\longrightarrow
\operatorname{erfc}
\longrightarrow
E\!\left[\frac{M_{\rm coll}}{M_0}\middle|M_0,z_0\right]
}.
$$

先给出最重要的结论：

$$
\boxed{
\text{是的，这个 EPS 解析式本质上给出一个条件期望的质量分数。}
}
$$

更严格地说，它不是一棵具体 merger tree 上确定的
$M_{\rm coll}^{\mathcal T}(z;f)/M_0$，而是

$$
\boxed{
E\!\left[
\frac{M_{\rm coll}^{\mathcal T}(z;f)}{M_0}
\,\middle|\,
\text{在 }z_0\text{ 形成质量 }M_0\text{ 的后代晕}
\right].
}
$$

Ludlow16 的记号把这个期望符号省略了。它用解析的平均 CMH 代表典型 CMH，
再用模拟重新校准密度关系的零点。

---

## 1. 先区分三种完全不同的“随机对象”

### 1.1 初始密度随机场

线性初始密度对比度写成

$$
\delta:\Omega\longrightarrow\mathcal S'(\mathbb R^3),
$$

其中

- $(\Omega,\mathcal F,P)$ 是底层概率空间；
- $\mathcal S'(\mathbb R^3)$ 是温和分布空间；
- $\delta(\omega)$ 是一次宇宙初始条件实现。

初始场的随机性来自概率测度 $P$，而不是来自空间坐标 $\mathbf x$。
空间点只是随机场的指标。

### 1.2 固定拉格朗日质量元的平滑轨迹

固定一个初始拉格朗日位置 $\mathbf q$，不断改变平滑尺度 $R$，得到

$$
R\longmapsto\delta_R(\mathbf q,\omega).
$$

对固定的 $\omega$，这是关于平滑尺度的一条轨迹。对固定的 $R$，
$\delta_R(\mathbf q)$ 是一个普通随机变量。

这里的“过程时间”不是宇宙时间，也不是红移，而是平滑尺度，后来会改用方差

$$
S\equiv\sigma^2(M)
$$

作为过程参数。

### 1.3 merger tree 上的 collapsed mass

给定一棵真实 merger tree $\mathcal T$，定义

$$
M_{\rm coll}^{\mathcal T}(z;f)
=
\sum_{i\in {\rm Anc}(\mathcal T,z)}
M_i(z)\,
\mathbf 1_{\{M_i(z)>fM_0\}}.
$$

这是 merger-tree 样本空间上的随机变量。不同初始场实现产生不同 merger tree，
也产生不同的 $M_{\rm coll}^{\mathcal T}$。

因此必须避免把下面三者混为一谈：

$$
\delta(\mathbf x),
\qquad
\delta(S),
\qquad
M_{\rm coll}^{\mathcal T}(z;f).
$$

第一个是空间随机场，第二个是尺度随机过程，第三个是由完整非线性演化和 halo finder
诱导出来的 merger-tree 泛函。

EPS 用第一和第二个对象近似第三个对象的条件平均。

---

## 2. 高斯随机场的概率测度

假设初始线性密度场是均值为零、统计均匀且各向同性的 Gaussian random field。
其 Fourier 协方差为

$$
E[\widetilde\delta(\mathbf k)
\widetilde\delta^*(\mathbf k')]
=(2\pi)^3\delta_D^{(3)}(\mathbf k-\mathbf k')P(k).
$$

$P(k)$ 是功率谱。它不是“另一个概率密度”，而是 Gaussian measure 的协方差算子
在 Fourier 空间中的表示。

对测试函数 $h\in\mathcal S(\mathbb R^3)$，定义线性观测

$$
\delta(h)=\langle\delta,h\rangle.
$$

高斯场的特征泛函为

$$
E[e^{i\delta(h)}]
=
\exp\!\left[-\frac12\langle h,\mathcal C_P h\rangle\right].
$$

在适当条件下，Bochner--Minlos 定理把这个特征泛函对应到
$\mathcal S'(\mathbb R^3)$ 上的高斯概率测度 $\mathcal G_P$：

$$
\mathcal G_P=\delta_*P.
$$

所以“抽取一组初始条件”可以写成

$$
\delta\sim\mathcal G_P.
$$

---

## 3. 平滑把随机场投影成一族随机变量

取窗口函数 $W_R$，定义

$$
\delta_R(\mathbf q)
=
\int d^3x\,
W_R(\mathbf q-\mathbf x)\delta(\mathbf x).
$$

Fourier 空间中是

$$
\delta_R(\mathbf q)
=
\int\frac{d^3k}{(2\pi)^3}
\widetilde W_R(k)\widetilde\delta(\mathbf k)e^{i\mathbf k\cdot\mathbf q}.
$$

因为这是 Gaussian random field 的线性泛函，任意有限组

$$
(\delta_{R_1}(\mathbf q),\ldots,\delta_{R_n}(\mathbf q))
$$

都联合高斯。因此

$$
\{\delta_R(\mathbf q):R>0\}
$$

是一个 Gaussian process。

它的协方差核为

$$
C(R_1,R_2)
=
E[\delta_{R_1}(\mathbf q)\delta_{R_2}(\mathbf q)]
=
\int_0^\infty\frac{k^2dk}{2\pi^2}
P(k)\widetilde W_{R_1}(k)\widetilde W_{R_2}(k).
$$

单点方差是

$$
S(R)=\sigma^2(R)
=
\int_0^\infty\frac{k^2dk}{2\pi^2}P(k)|\widetilde W_R(k)|^2.
$$

若用质量标记平滑尺度，则

$$
S(M)=\sigma^2(M).
$$

在通常 CDM 质量区间内，质量越小，纳入的高 $k$ 模式越多，因而

$$
M\downarrow
\quad\Longleftrightarrow\quad
S(M)\uparrow.
$$

---

## 4. 为什么 sharp-$k$ 窗口产生 Brownian motion

### 4.1 sharp-$k$ 平滑

取

$$
\widetilde W_K(k)=\mathbf 1_{\{k\le K\}}.
$$

于是

$$
\delta_K(\mathbf q)
=
\int_{|\mathbf k|\le K}
\frac{d^3k}{(2\pi)^3}
\widetilde\delta(\mathbf k)e^{i\mathbf k\cdot\mathbf q},
$$

并且

$$
S(K)
=
\int_0^K\frac{k^2dk}{2\pi^2}P(k).
$$

假设 $S(K)$ 严格递增，就可以反解 $K=K(S)$，定义

$$
B_S\equiv\delta_{K(S)}(\mathbf q).
$$

严格地说，sharp-$k$ 指示函数不是 Schwartz 测试函数，所以上式不应直接理解为
任意温和分布在普通测试函数上的配对。可以先在有限体积中离散 Fourier 模式，或把
$\delta_K(\mathbf q)$ 定义为均方意义下的随机积分；只要 $S(K)<\infty$，该随机积分就在
$L^2(\Omega,P)$ 中存在。随后利用下面的增量矩估计取得连续 modification，才得到
Wiener 路径空间中的随机元素。

### 4.2 协方差等于 $\min(S_1,S_2)$

若 $S_1\le S_2$，则 $K(S_1)\le K(S_2)$，故

$$
E[B_{S_1}B_{S_2}]
=
\int_0^{K(S_1)}\frac{k^2dk}{2\pi^2}P(k)
=S_1.
$$

所以一般地

$$
\boxed{E[B_{S_1}B_{S_2}]=\min(S_1,S_2).}
$$

### 4.3 增量独立

增量

$$
B_{S_2}-B_{S_1}
$$

只由 Fourier 壳层

$$
K(S_1)<|\mathbf k|\le K(S_2)
$$

中的模式贡献。不相交的 $S$ 区间对应不相交的 Fourier 壳层。

对于 Gaussian measure，协方差为零蕴含独立。因此

$$
B_{S_2}-B_{S_1}\sim N(0,S_2-S_1),
$$

而且不相交区间上的增量独立。

再加上 $B_0=0$ 和连续 modification 的存在，得到

$$
\boxed{\{B_S:S\ge0\}\text{ 是以 }S\text{ 为“时间”的标准 Brownian motion。}}
$$

这里 Brownian motion 描述的是：固定一个拉格朗日质量元，逐步提高质量分辨率时，
平滑密度对比度如何随机游走。它不是物质粒子在真实空间中的布朗运动。

---

## 5. Wiener measure 是怎样出现的

取路径空间

$$
\mathscr C=C_0([0,\infty),\mathbb R),
$$

先把“给定一个场实现、再在一个尺度上取值”写成二元映射。令
$u\in\mathcal S'(\mathbb R^3)$ 表示一个确定的场实现，$\mathsf S_K$ 表示
sharp-$k$ 平滑算子，定义

$$
\widehat\Phi_{\mathbf q}:
\mathcal S'(\mathbb R^3)\times[0,\infty)
\longrightarrow\mathbb R,
\qquad
\widehat\Phi_{\mathbf q}(u,S)
=
(\mathsf S_{K(S)}u)(\mathbf q).
$$

然后把第二个变量 $S$ 整体收进函数值中，得到与上式等价的路径值映射

$$
\Phi_{\mathbf q}:\mathcal S'(\mathbb R^3)\longrightarrow\mathscr C,
\qquad
[\Phi_{\mathbf q}(u)](S)
=
\widehat\Phi_{\mathbf q}(u,S).
$$

严格地说，选择连续 modification 后，这个式子先在一个满足
$\mathcal G_P(A_{\mathbf q})=1$ 的可测子集
$A_{\mathbf q}\subset\mathcal S'(\mathbb R^3)$ 上定义；再在零测补集上任意延拓，
便可把 $\Phi_{\mathbf q}$ 写成以上整个 $\mathcal S'$ 上的可测映射。这样的零测延拓
不会改变推前测度。

因此，$\widehat\Phi_{\mathbf q}$ 的确接收两个变量 $(u,S)$；
$\Phi_{\mathbf q}$ 只接收场实现 $u$，但返回整条函数
$S\mapsto\widehat\Phi_{\mathbf q}(u,S)$。两种写法的关系就是

$$
\boxed{
\widehat\Phi_{\mathbf q}(u,S)
=
e_S\!\left(\Phi_{\mathbf q}(u)\right),
}
$$

其中 $e_S$ 是下面定义的 evaluation map。

若初始随机场本身写成随机映射

$$
\delta:\Omega\longrightarrow\mathcal S'(\mathbb R^3),
$$

那么真正以基本样本 $\omega$ 为输入的路径值随机变量是

$$
\Psi_{\mathbf q}
=
\Phi_{\mathbf q}\circ\delta:
\Omega\longrightarrow\mathscr C,
$$

$$
[\Psi_{\mathbf q}(\omega)](S)
=
\widehat\Phi_{\mathbf q}(\delta(\omega),S).
$$

初始场高斯测度和底层概率测度满足

$$
\mathcal G_P=\delta_*P.
$$

所以路径分布可以等价地从两边推前：

$$
(\Psi_{\mathbf q})_*P
=
(\Phi_{\mathbf q})_*(\delta_*P)
=
(\Phi_{\mathbf q})_*\mathcal G_P.
$$

上一节已经证明其坐标过程具有 Brownian motion 的有限维分布和连续版本，因此

$$
\boxed{
(\Psi_{\mathbf q})_*P
=(\Phi_{\mathbf q})_*\mathcal G_P
=\mathbb W,
}
$$

其中 $\mathbb W$ 是路径空间上的 Wiener measure。

路径空间上的坐标映射是

$$
e_S:\mathscr C\to\mathbb R,
\qquad
e_S(\gamma)=\gamma(S).
$$

在 Wiener 空间

$$
(\mathscr C,\mathcal B(\mathscr C),\mathbb W)
$$

上，随机过程就是

$$
B_S=e_S.
$$

所以这里没有形式上的“对所有轨迹做普通积分”。概率是对 Wiener measure 积分：

$$
E_{\mathbb W}[F(B)]
=
\int_{\mathscr C}F(\gamma)\,d\mathbb W(\gamma).
$$

---

## 6. 从坍缩判据到首次穿越停时

### 6.1 红移给出吸收屏障

把线性扰动外推到今天，球对称坍缩阈值写成

$$
\delta_{\rm sc}(z)=\frac{\delta_c}{D(z)},
\qquad \delta_c\simeq1.686.
$$

在固定红移 $z$，这是 $B_S$ 路径空间中的常数屏障。

定义首次穿越停时

$$
\tau_{\delta}
=
\inf\{S\ge0:B_S\ge\delta\}.
$$

事件

$$
\tau_{\delta_{\rm sc}(z)}\in[S,S+dS]
$$

表示这个质量元第一次在方差尺度 $S$ 上进入一个于红移 $z$ 坍缩的区域。

“第一次”非常关键。如果只判断 $B_S>\delta$，同一个质量元可能同时被计入一个小晕
和包住它的大晕，这就是 cloud-in-cloud 问题。首次穿越把每条路径唯一分配给最先
跨越屏障的尺度。

### 6.2 后代晕是一个条件事件

考虑在观察红移 $z_0$、质量 $M_0$ 的后代晕。记

$$
S_0=S(M_0),
\qquad
\delta_0=\delta_{\rm sc}(z_0).
$$

在 excursion-set 语言中，这对应条件

$$
\tau_{\delta_0}=S_0.
$$

严格说，连续随机变量取精确值的事件通常概率为零。这里应理解为对
$\tau_{\delta_0}$ 的 regular conditional probability，或者先条件在
$[S_0,S_0+dS_0]$ 上再取 $dS_0\to0$。

这比简单条件 $B_{S_0}=\delta_0$ 更强，因为它还要求路径此前没有越过
$\delta_0$。

---

## 7. 强 Markov 性把条件问题重新变成 Brownian motion

取更早的红移 $z>z_0$。由于 $D(z)<D(z_0)$，有

$$
\delta_1\equiv\delta_{\rm sc}(z)>\delta_0.
$$

定义屏障差

$$
a\equiv\Delta\delta
=\delta_1-\delta_0>0.
$$

在后代首次穿越时刻 $S_0$ 以后，定义平移过程

$$
X_u=B_{S_0+u}-B_{S_0},
\qquad u\ge0.
$$

Brownian motion 在停时处具有强 Markov 性，因此在给定
$\tau_{\delta_0}=S_0$ 后，$X_u$ 在 regular conditional law 下仍是从零出发的
标准 Brownian motion，并且与 $S_0$ 以前的路径信息独立。

于是“质量元属于更早红移 $z$ 的某个祖先晕”变成：

> 从零出发的 Brownian path 是否在方差增量 $u$ 内首次到达高度 $a$？

定义

$$
T_a=\inf\{u\ge0:X_u\ge a\}.
$$

这就是一个标准的 Brownian first-passage problem。

---

## 8. 反射原理给出 $\operatorname{erfc}$

### 8.1 首次穿越概率

对从零出发的标准 Brownian motion，反射原理给出

$$
\mathbb W(T_a\le u)
=
\mathbb W\left(\sup_{0\le s\le u}X_s\ge a\right)
=2\,\mathbb W(X_u\ge a).
$$

由于

$$
X_u\sim N(0,u),
$$

所以

$$
\mathbb W(X_u\ge a)
=
\int_a^\infty
\frac{dx}{\sqrt{2\pi u}}
e^{-x^2/(2u)}
=
\frac12\operatorname{erfc}\!\left(\frac{a}{\sqrt{2u}}\right).
$$

乘上反射原理的因子 $2$，得到

$$
\boxed{
\mathbb W(T_a\le u)
=
\operatorname{erfc}\!\left(\frac{a}{\sqrt{2u}}\right).
}
$$

这解释了为什么最终答案不是

$$
\frac12\operatorname{erfc}(\cdots),
$$

而是完整的 $\operatorname{erfc}(\cdots)$。这个因子 $2$ 不是人为补丁；在
excursion-set 推导中，它来自首次穿越与反射原理，并自动解决 Press--Schechter
原始论证中的 cloud-in-cloud 归一化问题。

### 8.2 首次穿越时间的密度

对上式关于 $u$ 求导，得到 Lévy first-passage density：

$$
p_{T_a}(u)
=
\frac{a}{\sqrt{2\pi}}u^{-3/2}
\exp\!\left(-\frac{a^2}{2u}\right),
\qquad u>0.
$$

回到 EPS 记号，令

$$
u=S_1-S_0,
\qquad
a=\delta_1-\delta_0,
$$

便得到条件首次穿越密度

$$
\boxed{
f_{\rm FC}(S_1,\delta_1\mid S_0,\delta_0)
=
\frac{\delta_1-\delta_0}
{\sqrt{2\pi}(S_1-S_0)^{3/2}}
\exp\!\left[
-\frac{(\delta_1-\delta_0)^2}{2(S_1-S_0)}
\right].
}
$$

---

## 9. 从首次穿越分布到 progenitor mass fraction

### 9.1 方差与祖先质量方向相反

对早期祖先质量 $m\le M_0$，通常有

$$
S(m)\ge S(M_0)=S_0.
$$

阈值条件

$$
m\ge fM_0
$$

对应

$$
S_0\le S(m)\le S_f,
\qquad
S_f\equiv S(fM_0).
$$

因此允许的 Brownian “时间长度”为

$$
\Delta S_f=S_f-S_0
=\sigma^2(fM_0)-\sigma^2(M_0).
$$

### 9.2 一条 excursion-set 路径代表一个带标签的质量元

在后代晕的拉格朗日区域内均匀随机选一个质量元 $U$。令

$$
I_f(U,\mathcal T;z)
=
\mathbf 1_{\{U\text{ 在 }z\text{ 时属于质量 }m\ge fM_0\text{ 的祖先晕}\}}.
$$

对一棵固定 merger tree $\mathcal T$，在所有后代质量元上平均，有

$$
\frac1{M_0}\int_{{\cal H}_0}
I_f(U,\mathcal T;z)\,dM(U)
=
\frac{M_{\rm coll}^{\mathcal T}(z;f)}{M_0}.
$$

再对 merger-tree ensemble 取期望。由 Tonelli/Fubini 定理，

$$
E_{\mathcal T}\!\left[
\frac{M_{\rm coll}^{\mathcal T}(z;f)}{M_0}
\middle|M_0,z_0
\right]
=
P\!\left(
U\text{ 属于 }m\ge fM_0\text{ 的祖先}
\middle|M_0,z_0
\right).
$$

右边正是带标签质量元的 excursion-set 路径在 $\Delta S_f$ 以内穿越高屏障的条件概率：

$$
P(T_{\Delta\delta}\le\Delta S_f).
$$

这就是“概率”等于“平均质量分数”的数学原因。它依赖质量元的质量加权抽样，
不是把祖先晕按个数等权抽样。

### 9.3 与条件 progenitor mass function 的等价写法

若 $N(m,z\mid M_0,z_0)dm$ 表示平均祖先晕个数，则 EPS 给出

$$
\frac{dE[N]}{dm}
=
\frac{M_0}{m}
f_{\rm FC}(S(m),\delta_1\mid S_0,\delta_0)
\left|\frac{dS}{dm}\right|.
$$

乘以每个祖先贡献的质量 $m$，再除以 $M_0$，得到平均质量分数：

$$
E\!\left[\frac{M_{\rm coll}}{M_0}\middle|M_0,z_0\right]
=
\int_{fM_0}^{M_0}
\frac{m}{M_0}\frac{dE[N]}{dm}\,dm.
$$

变量代换 $m\mapsto S$ 后，

$$
=
\int_{S_0}^{S_f}
f_{\rm FC}(S_1,\delta_1\mid S_0,\delta_0)\,dS_1.
$$

所以这个公式求的是质量加权 conditional progenitor mass function 的累积，
不是祖先个数的期望。

---

## 10. 完整积分：得到 Ludlow16 的解析式

定义

$$
\Delta\delta
=
\delta_{\rm sc}(z)-\delta_{\rm sc}(z_0),
$$

$$
\Delta S_f
=
\sigma^2(fM_0)-\sigma^2(M_0).
$$

从首次穿越密度积分：

$$
E\!\left[
\frac{M_{\rm coll}^{\mathcal T}(z;f)}{M_0}
\middle|M_0,z_0
\right]
=
\int_0^{\Delta S_f}
\frac{\Delta\delta}{\sqrt{2\pi}}
u^{-3/2}
\exp\!\left[-\frac{(\Delta\delta)^2}{2u}\right]du.
$$

令

$$
y=\frac{\Delta\delta}{\sqrt{2u}},
$$

则

$$
du=-\frac{(\Delta\delta)^2}{y^3}dy,
$$

积分变为

$$
\frac{2}{\sqrt\pi}
\int_{\Delta\delta/\sqrt{2\Delta S_f}}^\infty
e^{-y^2}dy.
$$

按照 complementary error function 的定义

$$
\operatorname{erfc}(x)
=
\frac{2}{\sqrt\pi}\int_x^\infty e^{-y^2}dy,
$$

最终得到

$$
\boxed{
E\!\left[
\frac{M_{\rm coll}^{\mathcal T}(z;f)}{M_0}
\middle|M_0,z_0
\right]
=
\operatorname{erfc}\!\left[
\frac{\delta_{\rm sc}(z)-\delta_{\rm sc}(z_0)}
{\sqrt{2\left[\sigma^2(fM_0)-\sigma^2(M_0)\right]}}
\right].
}
$$

当条件集合把 $M_0,z_0$ 固定时，$M_0$ 是常数，因此也可写成

$$
\frac{E[M_{\rm coll}^{\mathcal T}(z;f)\mid M_0,z_0]}{M_0}
=\operatorname{erfc}(\cdots).
$$

但一般不能不加条件地写成

$$
E\!\left[\frac{M_{\rm coll}}{M_0}\right]
=
\frac{E[M_{\rm coll}]}{E[M_0]}.
$$

随机变量比值的期望通常不等于期望之比。这里只是因为模型先条件在固定的
$M_0,z_0$ 晕族群上，分母才可以移到期望外。

---

## 11. $\operatorname{erfc}$ 结果的三个极限检查

### 11.1 在观察红移处质量分数为一

当 $z=z_0$ 时，

$$
\Delta\delta=0,
$$

所以

$$
\operatorname{erfc}(0)=1.
$$

这表示在观察时刻，后代晕自身就是一个质量大于 $fM_0$ 的坍缩对象，故全部质量都计入。

### 11.2 向很早红移回溯时趋于零

当 $z\to\infty$ 时，通常

$$
\delta_{\rm sc}(z)\to\infty,
$$

从而

$$
\frac{M_{\rm coll}}{M_0}\to0.
$$

这表示极早时期几乎没有质量进入超过阈值 $fM_0$ 的坍缩祖先。

### 11.3 降低 $f$ 会增加 collapsed mass

当 $f$ 变小时，$fM_0$ 变小，通常 $S(fM_0)$ 变大，因此

$$
\Delta S_f\uparrow.
$$

$\operatorname{erfc}$ 的自变量减小，故

$$
E[M_{\rm coll}/M_0]\uparrow.
$$

这与“降低祖先质量门槛会统计更多坍缩结构”一致。

---

## 12. 它是均值、概率，还是中位数

### 12.1 在理想 EPS 模型中

三种说法等价：

$$
\boxed{
\begin{aligned}
&\text{带标签质量元在 }z\text{ 时属于 }m\ge fM_0\text{ 祖先的条件概率},\\
=\;&\text{条件 progenitor mass function 的累积质量分数},\\
=\;&E[M_{\rm coll}^{\mathcal T}/M_0\mid M_0,z_0].
\end{aligned}
}
$$

但它不等于：

- 一棵具体 merger tree 的确定 $M_{\rm coll}^{\mathcal T}/M_0$；
- 祖先晕数量 $N_{\rm prog}$ 的期望；
- 最大主祖先质量 $M_{\rm main}/M_0$ 的期望；
- $M_{\rm coll}/M_0$ 分布的中位数；
- 形成红移随机变量 $z_{-2}^{\mathcal T}$ 的期望。

### 12.2 “平均 CMH 的交点”不是“平均交点”

为了避免红移方向与时间方向相反造成歧义，先用宇宙时间定义个体晕的形成时刻：

$$
t_{-2}^{\mathcal T}
=
\inf\left\{t\le t_0:
M_{\rm coll}^{\mathcal T}(t;f)\ge M_{-2}^{\mathcal T}
\right\},
$$

再令

$$
z_{-2}^{\mathcal T}=z(t_{-2}^{\mathcal T}).
$$

解析模型则求

$$
E[M_{\rm coll}^{\mathcal T}(z;f)/M_0\mid M_0,z_0]
=
F_{-2}(c).
$$

由这个平均曲线得到的 $z_{-2}^{\rm EPS}$ 一般不满足

$$
z_{-2}^{\rm EPS}=E[z_{-2}^{\mathcal T}],
$$

也不必等于其中位数。因为“先取首次交点”和“先对曲线取平均再找交点”是两个
不可交换的非线性操作。

因此 Ludlow16 解析模型中的 $z_{-2}$ 更准确地说是一个有效的、族群层面的
characteristic collapse redshift。

---

## 13. 从平均 CMH 到浓度

给定剖面的内部质量比例

$$
F_{-2}(c)=\frac{M_{-2}}{M_0},
$$

Ludlow16 把解析平均 CMH 与内部特征质量匹配：

$$
F_{-2}(c)
=
\operatorname{erfc}\!\left[
\frac{\delta_{\rm sc}(z_{-2})-\delta_{\rm sc}(z_0)}
{\sqrt{2[S(fM_0)-S(M_0)]}}
\right].
$$

再联立密度记忆关系

$$
200c^3F_{-2}(c)
=
A_\rho
\frac{\rho_{\rm crit}(z_{-2})}{\rho_{\rm crit}(z_0)},
$$

求出有效的 $c(M_0,z_0)$ 和 $z_{-2}(M_0,z_0)$。

这一步的逻辑不是

$$
\text{一条 Brownian path}\longrightarrow\text{一个确定浓度},
$$

而是

$$
\boxed{
\text{Wiener path ensemble}
\longrightarrow
\text{平均 collapsed mass fraction}
\longrightarrow
\text{模拟校准的典型浓度关系}.
}
$$

若使用模拟中每个个体晕的真实 merger tree，则可以直接构造
$M_{\rm coll}^{\mathcal T}$，不必用 EPS 的平均 erfc 曲线。Ludlow16 对这两种情况使用
不同的密度归一化常数：真实 merger tree 约为 $A_\rho\simeq400$，解析 EPS 模型则
采用约 $A_\rho\simeq650$，用来补偿球对称坍缩和解析 CMH 的系统误差。

---

## 14. 精确 Brownian 推导与 Ludlow16 实际实现之间的差别

这一节非常重要。上面的 Brownian first-crossing 推导是严格的，但严格性依赖一组
理想化假设。

### 14.1 sharp-$k$ 才给出独立增量

精确的

$$
E[B_{S_1}B_{S_2}]=\min(S_1,S_2)
$$

来自 sharp-$k$ 窗口。

对 real-space top-hat 窗口，两个平滑尺度使用的 Fourier 模式权重彼此重叠，过程满足

$$
E[\delta(S_1)\delta(S_2)]
\ne\min(S_1,S_2),
$$

增量相关，过程一般不是 Markov，也不是标准 Brownian motion。此时简单的反射原理
和闭式 $\operatorname{erfc}$ 不再严格成立。

Bond et al. (1991) 明确指出：sharp-$k$ 情形是 Brownian random walk，首次穿越可解析；
Gaussian 或 real-space top-hat 等一般滤波需要处理相关步长，通常没有同样简单的闭式解。

### 14.2 Ludlow16 使用的是常见的混合近似

Ludlow16 实际从线性功率谱用 real-space spherical top-hat 计算 $\sigma(M)$，但仍采用
标准 Markov EPS 的 $\operatorname{erfc}$ 形式。因此应理解为：

$$
\boxed{
\text{Brownian/sharp-}k\text{ 推导提供公式结构，}
\quad
\text{top-hat }\sigma(M)\text{ 与模拟校准提供现象学修正。}
}
$$

所以不能声称 Ludlow16 的实际 top-hat 实现严格等于 Wiener measure 下的首次穿越概率；
它是由理想 Markov excursion set 启发、再经过数值校准的近似。

### 14.3 常数屏障也是近似

闭式结果还假设：

- 初始场为 Gaussian；
- 坍缩可由常数球对称屏障描述；
- 过程是 sharp-$k$ Markov walk；
- halo 质量能与首次穿越尺度一一对应；
- 不显式处理潮汐、椭球坍缩、环境和 assembly bias。

移动屏障、非高斯初始条件或相关步长都会改变 first-crossing law。

Ludlow16 在展示 CMH 形状时还尝试过低于标准 $1.686$ 的有效坍缩阈值，以补偿球对称
坍缩误差；这再次说明实际模型包含模拟校准，而不是纯数学定理直接给出全部数值。

---

## 15. 它为什么能够包含“很多块祖先晕”

### 15.1 首次穿越到宿主质量，是 excursion-set 的质量分配公设

固定一个拉格朗日质量元 $\mathbf q$ 和红移 $z$，令

$$
\mathcal E_z(\mathbf q)
=
\left\{M:\delta_M(\mathbf q)\geq\delta_{\rm sc}(z)\right\}.
$$

由于 $S(M)$ 随 $M$ 减小，首次穿越时间

$$
\tau_z(\mathbf q)
=
\inf\left\{S:\delta_S(\mathbf q)\geq\delta_{\rm sc}(z)\right\}
$$

对应满足坍缩判据的最大平滑质量：

$$
M_{\rm ES}(\mathbf q,z)
=M[\tau_z(\mathbf q)]
=\sup\mathcal E_z(\mathbf q).
$$

excursion-set 随后**规定**把该质量元分配给质量为 $M_{\rm ES}$ 的宿主晕。因此

$$
M_{\rm host}^{\rm(real)}(\mathbf q,z)
\simeq M_{\rm ES}(\mathbf q,z)
$$

是模型的统计对应，而不是由 Brownian motion 严格推出的动力学恒等式。首次穿越只给出
候选宿主的质量尺度，不能单独确定真实晕的中心、边界、形状或成员集合。Wiener 测度使
这套分配规则的概率可解析计算，但不会使它成为关于真实非线性晕的严格定理。

若先在模拟中识别了真实后代晕，严格条件应是“$\mathbf q$ 属于该晕且该晕质量为
$M_0$”；EPS 则用 $\{\tau_{\delta_0}=S_0\}$ 代替这个条件。两者不是事件集合上的恒等，
而是解析模型对真实选晕条件的近似编码。

### 15.2 一条路径由 $(\omega,\mathbf q)$ 共同确定

完整的路径值映射应写成

$$
\Gamma:\Omega\times\mathbb R^3_{\rm L}\longrightarrow\mathscr C,
\qquad
[\Gamma(\omega,\mathbf q)](S)=\delta_S(\mathbf q;\omega).
$$

所以，“一条路径代表一个宇宙实现 $\omega$”和“一条路径代表一个带标记质量元
$\mathbf q$”并不字面等价：

- 固定 $\mathbf q$、改变 $\omega$：得到该位置在不同宇宙实现中的路径系综；在
  sharp-$k$ 情形下，其推前测度是 Wiener 测度。
- 固定 $\omega$、改变 $\mathbf q$：得到同一宇宙实现中不同质量元的路径族；这些路径
  通常彼此相关。

只有再假设统计齐次性与遍历性，固定 $\mathbf q_0$ 的系综平均才等于单个无限大实现中
的空间（初始质量）平均。对合适的单路径泛函 $F$，

$$
\mathbb E_{\mathbb P}
\!\left[F(\Gamma(\omega,\mathbf q_0))\right]
=
\lim_{V\to\infty}\frac{1}{|V|}
\int_V F(\Gamma(\omega,\mathbf q))\,d^3q
\quad\text{a.s.}
$$

这里相等的是两种**平均**，不是两个标签本身。单点 Wiener 测度也不包含不同质量元
如何组合成同一个真实晕的联合空间信息。

### 15.3 单路径平均为什么能给出多祖先总质量

一条 excursion-set path 不是一棵完整 merger tree。它只描述随机选中的一个后代质量元
在不同质量尺度上的嵌套环境。

但对所有质量元取平均时：

- 落在祖先 $A$ 中的质量元贡献 $M_A/M_0$；
- 落在祖先 $B$ 中的质量元贡献 $M_B/M_0$；
- 依此类推。

因此概率自动产生

$$
\frac{M_A+M_B+\cdots}{M_0},
$$

即所有超过阈值祖先的总质量分数。这就是为什么单质量元 random walk 能给出
多祖先总质量的期望，却不能告诉我们祖先的具体数目和组合方式。

例如它可以预测

$$
E[M_{\rm coll}/M_0]=0.4,
$$

但无法仅凭这个数区分：

$$
0.4=0.4,
\qquad
0.4=0.2+0.2,
\qquad
0.4=0.15+0.15+0.10.
$$

这些 merger configurations 对 CMH 质量分数相同，但动力学状态可以完全不同。

---

## 16. 与“同一个晕跨红移自洽”的关系

EPS 公式在每个 $(M_0,z_0)$ 上给出一个 regular conditional ensemble：

$$
\mu_{M_0,z_0}
=
\mathbb W(\,\cdot\mid\tau_{\delta_0}=S_0).
$$

改变观察红移或终点质量，就改变条件测度：

$$
\mu_{M_0,z_0}
\longrightarrow
\mu_{M_1,z_1}.
$$

因此分别计算

$$
\bar c(M_0,z_0),
\qquad
\bar c(M_1,z_1)
$$

并不会自动生成同一条个体 halo path 的联合分布。它只保证每个端点条件下的边缘统计
关系正确或近似正确。

若要求同一个晕跨快照严格自洽，需要额外构造：

- 一棵共享同一后代--祖先关系的 Monte Carlo EPS merger tree；或
- 一棵由 $N$-body 粒子追踪得到的真实 merger tree。

然后在同一棵树上定义所有 $M_{\rm coll}^{\mathcal T}(z;f)$。单独的 erfc 公式不包含
这种跨红移的联合耦合信息。

---

## 17. 一张总表

| 层次 | 数学对象 | 概率测度/分布 | 物理含义 |
|---|---|---|---|
| 初始条件 | $\delta\in\mathcal S'(\mathbb R^3)$ | Gaussian measure $\mathcal G_P$ | 一次初始密度场实现 |
| 尺度过程 | $S\mapsto\delta(S)$ | sharp-$k$ 时为 Wiener measure $\mathbb W$ | 固定质量元随分辨率变化的平滑密度 |
| 坍缩判据 | $\tau_\delta$ | first-passage law | 第一次跨越坍缩屏障的尺度 |
| 条件后代 | $\tau_{\delta_0}=S_0$ | regular conditional law | EPS 中质量元被分配到 $(M_0,z_0)$ 后代晕；是真实选晕条件的近似编码 |
| 祖先分布 | $f_{\rm FC}(S_1\mid S_0)$ | 条件首次穿越密度 | 质量加权 progenitor distribution |
| CMH | $M_{\rm coll}^{\mathcal T}/M_0$ | merger-tree ensemble 上的随机变量 | 所有大于 $fM_0$ 祖先的总质量分数 |
| EPS 解析式 | $\operatorname{erfc}(\cdots)$ | 条件 first-crossing 概率 | $E[M_{\rm coll}/M_0\mid M_0,z_0]$ |
| Ludlow 浓度 | $\bar c(M_0,z_0)$ | 解析平均加模拟校准 | 弛豫晕族群的典型浓度 |

---

## 18. 最终物理与测度论图像

最简洁的测度论表达是：初始高斯场测度 $\mathcal G_P$ 经由 sharp-$k$ 平滑路径映射
$\Phi_{\mathbf q}$ 推前为 Wiener measure，

$$
(\Phi_{\mathbf q})_*\mathcal G_P=\mathbb W.
$$

在 Wiener 路径空间上，halo formation 被近似为常数屏障的首次穿越。给定后代首次穿越
$(S_0,\delta_0)$ 后，强 Markov 性允许把剩余路径重新视为从零开始的 Brownian motion。
更早屏障与后代屏障之差为

$$
\Delta\delta
=\delta_{\rm sc}(z)-\delta_{\rm sc}(z_0),
$$

祖先阈值给出的最大允许方差增量为

$$
\Delta S_f
=\sigma^2(fM_0)-\sigma^2(M_0).
$$

反射原理于是给出

$$
\mathbb W(T_{\Delta\delta}\le\Delta S_f)
=
\operatorname{erfc}\!\left(
\frac{\Delta\delta}{\sqrt{2\Delta S_f}}
\right).
$$

因为 excursion-set path 对应随机抽取的后代质量元，这个路径概率等于该质量元落入
超过阈值祖先的概率；再由 Tonelli/Fubini，它等于整个后代晕中这类质量元所占比例的
条件期望：

$$
\boxed{
\mathbb W(T_{\Delta\delta}\le\Delta S_f)
=
E\!\left[
\frac{M_{\rm coll}^{\mathcal T}(z;f)}{M_0}
\middle|M_0,z_0
\right].
}
$$

这就是 Ludlow16 的 $\operatorname{erfc}$ 解析式从高斯随机场、尺度随机过程、
Brownian motion 和 Wiener measure 出发的完整含义。

---

## 19. 参考文献与适用范围

1. J. R. Bond, S. Cole, G. Efstathiou & N. Kaiser,
   [Excursion set mass functions for hierarchical Gaussian fluctuations](https://articles.adsabs.harvard.edu/pdf/1991ApJ...379..440B),
   *ApJ* **379**, 440 (1991)。该文建立 excursion-set 首次穿越框架，并明确说明
   sharp-$k$ 滤波产生 Brownian random walk，而一般滤波导致相关步长。
2. C. Lacey & S. Cole,
   [Merger rates in hierarchical models of galaxy formation](https://adsabs.harvard.edu/pdf/1993MNRAS.262..627L),
   *MNRAS* **262**, 627 (1993)。该文系统建立后代条件下的 progenitor mass function
   和 EPS merger history。
3. A. D. Ludlow et al.,
   [The mass--concentration--redshift relation of cold and warm dark matter haloes](https://academic.oup.com/mnras/article/460/2/1214/2608988),
   *MNRAS* **460**, 1214 (2016)。其 equation (3) 使用上述累积条件质量分数，
   equations (6)--(7) 再把 characteristic collapse redshift 与 halo concentration 连接。

本文的严格 Brownian/Wiener 推导适用于 Gaussian 初始场、sharp-$k$ Markov walk 和
常数球对称坍缩屏障。Ludlow16 的实际数值模型使用 real-space top-hat 计算
$\sigma(M)$，并通过模拟校准补偿这些理想化假设，因此其最终浓度关系应理解为经过验证的
统计模型，而不是 Wiener measure 单独推出的无参数定理。
