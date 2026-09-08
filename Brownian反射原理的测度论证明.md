# Brownian 反射原理的测度论证明

本文补充《Ludlow16-CMH解析式的测度论推导》第 8 节使用的反射原理，目标是严格说明

$$
P(T_a\le t)
=
2P(B_t\ge a)
=
\operatorname{erfc}\!\left(\frac{a}{\sqrt{2t}}\right),
$$

以及因子 $2$ 为什么来自路径空间上的保测反射，而不是人为加入的修正。

---

## 1. 概率空间、坐标过程与首次到达时间

取 canonical Wiener space

$$
(\mathscr C,\mathcal B(\mathscr C),\mathbb W),
\qquad
\mathscr C=C_0([0,\infty),\mathbb R).
$$

路径空间的元素记为

$$
\gamma:[0,\infty)\to\mathbb R,
\qquad \gamma(0)=0.
$$

坐标过程定义为

$$
B_t(\gamma)=\gamma(t).
$$

在 Wiener measure $\mathbb W$ 下，$B_t$ 是标准 Brownian motion。固定屏障 $a>0$，
定义首次到达时间

$$
T_a(\gamma)
=
\inf\{s\ge0:\gamma(s)=a\}.
$$

由于路径连续，

$$
\{T_a\le t\}
=
\left\{\sup_{0\le s\le t}B_s\ge a\right\}.
$$

$T_a$ 是关于自然滤过

$$
\mathcal F_t=\sigma(B_s:0\le s\le t)
$$

的停时，因为事件 $\{T_a\le t\}$ 只依赖截至 $t$ 的路径。

---

## 2. 在首次碰撞后反射路径

对满足 $T_a(\gamma)<\infty$ 的路径，定义反射映射

$$
(R_a\gamma)(s)
=
\begin{cases}
\gamma(s),&0\le s\le T_a(\gamma),\\[4pt]
2a-\gamma(s),&s>T_a(\gamma).
\end{cases}
$$

因为

$$
\gamma(T_a)=a,
$$

所以反射后的左右极限都等于 $a$，$R_a\gamma$ 仍是连续路径。

反射两次恢复原路径：

$$
R_a(R_a\gamma)=\gamma.
$$

因此

$$
R_a^{-1}=R_a.
$$

它是路径集合上的一一对应，而不仅是终点数值的变换。

---

## 3. 反射为什么保持 Wiener measure

在停时 $T_a$ 后定义增量过程

$$
X_u
=
B_{T_a+u}-B_{T_a},
\qquad u\ge0.
$$

Brownian motion 的强 Markov 性说明，在 $\{T_a<\infty\}$ 上，条件于
$\mathcal F_{T_a}$ 后：

1. $X$ 是从零出发的标准 Brownian motion；
2. $X$ 与碰撞以前的信息 $\mathcal F_{T_a}$ 独立。

而标准 Brownian motion 关于零对称：

$$
X\overset{d}{=}-X.
$$

反射后的路径在 $T_a$ 之后为

$$
(R_aB)_{T_a+u}
=
2a-B_{T_a+u}
=
a-X_u.
$$

原路径在该处为

$$
B_{T_a+u}=a+X_u.
$$

所以 $R_a$ 在停时以后所做的事情恰好是

$$
X\longmapsto-X.
$$

由于 $X$ 与 $-X$ 具有相同 Wiener 分布，反射不会改变停时后的条件路径律。

更严格地，若 $G$ 是碰撞前路径的有界 $\mathcal F_{T_a}$-可测泛函，$H$ 是碰撞后
增量路径的有界 Borel 泛函，则

$$
E_{\mathbb W}\!\left[
G\,H(X)
\right]
=
E_{\mathbb W}\!\left[
G\,H(-X)
\right].
$$

证明为

$$
\begin{aligned}
E_{\mathbb W}[G H(X)]
&=
E_{\mathbb W}\!\left[
G\,E_{\mathbb W}[H(X)\mid\mathcal F_{T_a}]
\right]\\
&=
E_{\mathbb W}[G]\,E_{\mathbb W}[H(B)]\\
&=
E_{\mathbb W}[G]\,E_{\mathbb W}[H(-B)]\\
&=
E_{\mathbb W}[G H(-X)].
\end{aligned}
$$

通过单调类定理从乘积型泛函推广到路径空间上的全部相关 Borel 事件，得到

$$
\boxed{
\mathbb W(R_a^{-1}A)=\mathbb W(A)
}
$$

对反射所涉及的路径事件成立。这就是“反射保持 Wiener measure”的严格含义。

---

## 4. 两个路径事件之间的双射

固定 $b<a$，定义

$$
A_{a,b}
=
\{T_a\le t,\ B_t\le b\},
$$

$$
C_{a,b}
=
\{B_t\ge2a-b\}.
$$

若 $\gamma\in A_{a,b}$，则反射后终点满足

$$
(R_a\gamma)(t)
=
2a-\gamma(t)
\ge2a-b.
$$

所以

$$
R_a(A_{a,b})\subseteq C_{a,b}.
$$

反过来，若 $\gamma\in C_{a,b}$，则

$$
\gamma(t)\ge2a-b>a.
$$

因为 $\gamma(0)=0<a$ 且路径连续，它在 $t$ 以前必定经过 $a$，所以 $T_a\le t$。
反射后的终点满足

$$
(R_a\gamma)(t)
=
2a-\gamma(t)
\le b.
$$

因此

$$
R_a(C_{a,b})\subseteq A_{a,b}.
$$

结合 $R_a^{-1}=R_a$，得到双射

$$
\boxed{
R_a:A_{a,b}\longleftrightarrow C_{a,b}.
}
$$

再利用保测性：

$$
\boxed{
P(T_a\le t,B_t\le b)
=
P(B_t\ge2a-b),
\qquad b<a.
}
$$

这就是反射原理的一般形式。

---

## 5. 从一般形式得到因子 $2$

把已经碰过屏障的路径按终点分成两类：

$$
\{T_a\le t\}
=
\{B_t\ge a\}
\mathbin{\dot\cup}
\{T_a\le t,B_t<a\}.
$$

第一类路径在终点仍位于屏障之上。第二类路径曾经碰到屏障，但终点已经落回屏障以下。

在上一节令 $b\uparrow a$。由概率测度的下连续性，

$$
P(T_a\le t,B_t<a)
=
P(B_t>a).
$$

因为 $B_t\sim N(0,t)$ 具有连续密度，

$$
P(B_t=a)=0,
$$

所以

$$
P(B_t>a)=P(B_t\ge a).
$$

于是

$$
\begin{aligned}
P(T_a\le t)
&=P(B_t\ge a)+P(T_a\le t,B_t<a)\\
&=P(B_t\ge a)+P(B_t>a)\\
&=2P(B_t\ge a).
\end{aligned}
$$

因此

$$
\boxed{
P(T_a\le t)=2P(B_t\ge a).
}
$$

---

## 6. 从高斯尾概率得到 $\operatorname{erfc}$

由于

$$
B_t\sim N(0,t),
$$

有

$$
P(B_t\ge a)
=
\int_a^\infty
\frac{1}{\sqrt{2\pi t}}
\exp\!\left(-\frac{x^2}{2t}\right)dx.
$$

令

$$
y=\frac{x}{\sqrt{2t}},
$$

得到

$$
P(B_t\ge a)
=
\frac12\operatorname{erfc}\!\left(
\frac{a}{\sqrt{2t}}
\right).
$$

代入反射原理：

$$
\boxed{
P(T_a\le t)
=
\operatorname{erfc}\!\left(
\frac{a}{\sqrt{2t}}
\right).
}
$$

对 $t$ 求导还可得到首次到达时间密度

$$
\boxed{
p_{T_a}(t)
=
\frac{a}{\sqrt{2\pi}t^{3/2}}
\exp\!\left(-\frac{a^2}{2t}\right),
\qquad t>0.
}
$$

---

## 7. 在 EPS/Ludlow16 中的变量替换

给定后代首次穿越点

$$
(S_0,\delta_0)
=
(S(M_0),\delta_{\rm sc}(z_0)),
$$

强 Markov 性允许在该点重新开始。令

$$
a
=
\Delta\delta
=
\delta_{\rm sc}(z)-\delta_{\rm sc}(z_0),
$$

$$
t
=
\Delta S_f
=
\sigma^2(fM_0)-\sigma^2(M_0).
$$

于是

$$
P(T_{\Delta\delta}\le\Delta S_f)
=
\operatorname{erfc}\!\left[
\frac{\delta_{\rm sc}(z)-\delta_{\rm sc}(z_0)}
{\sqrt{2[\sigma^2(fM_0)-\sigma^2(M_0)]}}
\right].
$$

这给出 EPS 中随机质量元落入质量 $m\ge fM_0$ 祖先的条件概率，进一步通过质量元平均
与 merger-tree ensemble 平均，成为 Ludlow16 使用的平均 collapsed mass fraction。

---

## 8. 反射原理依赖哪些假设

上述证明使用了：

1. 路径连续，从终点越过 $a$ 可以推出此前必经 $a$；
2. $T_a$ 是停时；
3. Brownian motion 在停时处具有强 Markov 性；
4. Brownian 增量关于零对称，$B\overset d=-B$；
5. $B_t$ 没有原子，因此 $P(B_t=a)=0$。

在 excursion-set 理论中，这些性质严格对应 Gaussian 初始场与 sharp-$k$ 滤波产生的
Markov Brownian walk。real-space top-hat 滤波产生相关步长时，不能未经修正地直接使用
这套反射论证。

---

## 9. 一句话总结

反射映射把“在 $t$ 前碰过 $a$、但终点落回 $a$ 以下”的路径，与“终点位于 $a$ 以上”
的路径逐一配对；强 Markov 性和 Brownian 对称性保证每对路径集合具有相同 Wiener
测度。因此首次越障概率等于终点高斯尾概率的两倍：

$$
\boxed{
\mathbb W(T_a\le t)
=2\mathbb W(B_t\ge a)
=\operatorname{erfc}\!\left(\frac{a}{\sqrt{2t}}\right).
}
$$
