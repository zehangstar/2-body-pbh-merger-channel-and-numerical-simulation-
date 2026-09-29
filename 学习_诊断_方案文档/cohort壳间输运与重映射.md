# cohort 壳间输运 / 重映射方案

> 状态：notebook 级实验方案。  
> 当前使用位置：`M_010_15j_then_a.ipynb`、`M_010_15a_then_j.ipynb`。  
> 当前未写入：`binary_single_montecarlo.py`、`binary_single.py`、`halostructure.py` 等主体程序。  
> 解释边界：这是为方案 B 补充的项目闭合条件，不是目标论文已经公开的壳间输运算法。

## 1. 为什么需要壳间输运

方案 B 把每个壳中的存活双星 cohort 从一个全局时间边界继承到下一个边界。原始实现默认把各壳质量正增量

$$
\Delta M_i^{(k)}=M_i^{(k)}-M_i^{(k-1)}
$$

都解释为进入第 $i$ 壳的新吸积物质，并要求 $\Delta M_i^{(k)}\ge 0$。这对低质量 halo 的既有计算可以工作，但不构成一般的径向输运模型。

当 $M_0=10^{15}M_\odot$、`prada12_hmf`、从 $z=12$ 开始、全局环境步长为 200 Myr 且划分 10 壳时，当前壳结构序列出现：

- 69 个时间边界、68 个演化区间；
- 690 个“边界—壳”质量单元中有 39 个负增量；
- 负增量全部出现在第 2 壳；
- 最负的单步增量约为 $-3.93\times10^3M_\odot$，约为该壳质量的 $-0.4\%$；
- 同期 halo 总质量在所有时间步都增长。

因此，负增量不是 halo 整体失去质量，而是随 NFW 结构、浓度和壳边界变化产生的壳间质量重分配。若直接报错，方案 B 不能继续；若把负值截断为零或取绝对值，则会凭空创造质量并破坏 cohort 人口账本。本方案显式加入壳间输运，保留并暴露这种重分配。

## 2. 本方案解决什么、不解决什么

本方案解决的是数值人口闭合问题：

1. 将负壳质量增量解释为旧物质离开 donor 壳；
2. 将这部分旧物质分配给有正质量增量的 receiver 壳；
3. receiver 尚未由旧物质填充的正增量才解释为新的 halo 吸积；
4. 将存活 tracer 的轨道状态和人口权重一起重映射；
5. 在人口守恒式中显式加入输运流入和流出。

本方案不声称求解真实的径向相空间动力学。它没有为每个双星积分半径、径向速度、动力摩擦或扩散方程，也没有从论文中恢复出作者未公开的输运实现。

## 3. 符号与状态

在时间边界 $k$，第 $i$ 壳的环境质量记为

$$
M_i^{(k)}.
$$

第 $i$ 壳的存活 tracer 集合记为

$$
\mathcal P_i^{(k)}=
\left\{a_n,e_n,w_n,b_n\right\},
$$

其中：

- $a_n$：半长轴；
- $e_n$：偏心率；
- $w_n$：该 tracer 代表的物理双星人口权重；
- $b_n$：其首次进入 halo 的时间边界 `entry_boundary`。

壳间迁移时保留这四项，不重新抽样轨道，也不重置 cohort 身份。

## 4. 第一步：确定性的质量流闭合

相邻边界的壳质量变化为

$$
\Delta M_i^{(k)}=M_i^{(k)}-M_i^{(k-1)}.
$$

定义 donor 需求和 receiver 容量：

$$
D_i^{(k)}=\max\!\left(-\Delta M_i^{(k)},0\right),
\qquad
G_i^{(k)}=\max\!\left(\Delta M_i^{(k)},0\right).
$$

在 halo 总质量增长时，必有

$$
\sum_iG_i^{(k)}\ge \sum_iD_i^{(k)}.
$$

输运矩阵

$$
T_{ij}^{(k)}\ge0
$$

表示从旧边界的 donor 壳 $i$ 转移到新边界 receiver 壳 $j$ 的环境质量。实现名为

```text
nearest_shell_minimum_distance_v1
```

其精确规则是：

1. 列出所有 `donor -> receiver` 组合；
2. 按壳编号距离 $|i-j|$ 从小到大排序；
3. 距离相同时按 receiver 编号、再按 donor 编号确定顺序；
4. 每次输送 donor 剩余需求与 receiver 剩余容量中的较小者；
5. 直到全部 donor 需求被吸收。

因此它是确定性的“最近壳优先”分配。在目前 $M_0=10^{15}M_\odot$ 的轨迹中，每个异常边界只有第 2 壳是 donor，所以不存在多个 donor 竞争造成的歧义。名称中的 `minimum_distance` 表示最近壳优先的工程闭合，不应解释成已证明的一般全局最优输运求解器。

输运约束为

$$
\sum_jT_{ij}^{(k)}=D_i^{(k)},
\qquad
\sum_iT_{ij}^{(k)}\le G_j^{(k)}.
$$

receiver 壳中剩余的新吸积质量定义为

$$
A_j^{(k)}=G_j^{(k)}-\sum_iT_{ij}^{(k)}\ge0.
$$

于是每个壳都满足确定性的质量重建式

$$
M_i^{(k)}=
M_i^{(k-1)}
-\sum_jT_{ij}^{(k)}
+\sum_jT_{ji}^{(k)}
+A_i^{(k)}.
$$

初始边界没有上一时刻，约定

$$
A_i^{(0)}=M_i^{(0)},\qquad T_{ij}^{(0)}=0.
$$

当前 10 壳测试中，39 个负增量单元对应的计划输运总质量约为

$$
\sum_{kij}T_{ij}^{(k)}=1.0437139797\times10^5M_\odot.
$$

## 5. 第二步：存活 tracer 的随机重映射

质量流 $T_{ij}^{(k)}$ 是确定的，但现有 cohort 只保存加权 tracer，没有保存每个双星的连续径向位置。因此必须把环境质量流转化为 tracer 迁移概率。

对 donor 壳 $i$ 中的每个存活 tracer，定义

$$
p_{i\to j}^{(k)}=
\frac{T_{ij}^{(k)}}{M_i^{(k-1)}},
$$

以及留在原壳的概率

$$
p_{i\to i}^{(k)}=
1-\sum_{j\ne i}p_{i\to j}^{(k)}.
$$

使用当前作业的 NumPy 随机数发生器对每个 tracer 抽取一个且仅一个目标壳。若 tracer 从 $i$ 迁移到 $j$，则完整保留

$$
(a,e,w,b)_{i}\longrightarrow(a,e,w,b)_{j}.
$$

此处有三个重要性质：

- 不复制 tracer，因此不会由于输运重复计算同一物理人口；
- 不拆分单个 tracer 的权重，因此不会造成 tracer 数随时间指数增长；
- 不重新采样 $(a,e)$，因此 cohort 已经积累的轨道演化不会丢失。

由于这是 Monte Carlo 实现，实际迁移的存活权重

$$
W_{ij}^{(k)}=
\sum_{n\in i\to j}w_n
$$

会围绕质量流所暗示的期望值波动。程序的人口账本使用实际实现的 $W_{ij}^{(k)}$，而不是用期望值代替它。样本量增大时，这一重映射噪声应下降；它必须作为单独的数值误差来源进行收敛检查。

## 6. 第三步：只对真正的新吸积质量生成 cohort

定义单位环境质量对应的原初双星系统数

$$
\eta_b=
\frac{f_{\rm PBH}f_{\rm binary}}{m_1+m_2}.
$$

边界 $k$、壳 $i$ 的新供给物理权重为

$$
Q_i^{(k)}=\eta_b A_i^{(k)}.
$$

若 $Q_i^{(k)}>0$，则：

- 初始边界抽取 `N1` 个 tracer；
- 后续边界抽取 `N2` 个 tracer；
- 每个新 tracer 的权重为 $Q_i^{(k)}/N$。

若 $A_i^{(k)}=0$，该壳在这一边界不抽取 `N2`。所以 notebook 中给出的

$$
N_{\rm upper}=N_{\rm shell}
\left(N_1+N_{\rm boundary}N_2\right)
$$

只是保守上限，最终以 `drawn_per_boundary.sum()` 或汇总文件中的 `samples_drawn` 为准。

## 7. `gw_aged` 入口与单边界更新顺序

本方案保留原方案 B 的 `gw_aged` 入口。新抽取的原初双星先从设定形成时间通过孤立 Peters GW 演化到入晕边界：

1. 若在入晕前已经并合，其权重计入 `outside_merged_weight_per_boundary`，不进入 halo；
2. 未并合者以入晕时的 $(a,e)$ 状态加入相应壳；
3. 已在 halo 中的旧 cohort 不重复做入晕前老化。

时间边界 $k$ 的完整顺序是：

```text
读取边界 k 的壳环境
    -> 对上一步存活人口执行壳间重映射
    -> 按 A_i(k) 抽取并老化新 cohort
    -> 合并旧人口与新人口
    -> 记录当前边界存活权重
    -> 在边界 k 的冻结环境中推进到 k+1
    -> 永久移除晕内并合者
    -> 写入检查点
```

最后一个 $z=0$ 边界仍登记输运和新供给，但没有额外的后续演化时长，与原方案 B 的边界语义一致。

## 8. 加入输运后的逐壳人口守恒

对第 $i$ 壳，在边界 $k$ 定义：

- $Q_i$：新吸积供给权重；
- $W_{ji}$：从其他壳迁入的实际存活权重；
- $W_{ij}$：迁往其他壳的实际存活权重；
- $O_i$：入晕前已经并合的权重；
- $I_i$：入晕后已经并合的权重；
- $S_i$：当前仍存活的权重。

程序检查

$$
\mathcal R_i^{(k)}=
\sum_{q\le k}
\left[
Q_i^{(q)}+
\sum_jW_{ji}^{(q)}-
\sum_jW_{ij}^{(q)}-
O_i^{(q)}
\right]
-\sum_{q<k}I_i^{(q)}
-S_i^{(k)}.
$$

理想情况下 $\mathcal R_i^{(k)}=0$，实际只保留浮点舍入误差。全 halo 求和时，所有壳间迁移项严格成对抵消：

$$
\sum_i\sum_jW_{ji}
-\sum_i\sum_jW_{ij}=0.
$$

因此全局守恒仍是

$$
\text{累计新供给}
=\text{晕外并合}
+\text{晕内并合}
+\text{当前存活}.
$$

输出同时报告：

- `shell_balance_max_abs`；
- `global_balance_max_abs`；
- `shell_balance_max_relative`；
- `global_balance_max_relative`。

绝对残差应结合总人口权重理解；在 $10^{15}M_\odot$ halo 中权重很大，约 $10^{-3}$ 的绝对舍入差仍可能对应约 $10^{-16}$ 的相对误差。

## 9. 保存的输运与事件对象

实验历史对象在原 `CohortShellMergerHistory` 基础上增加：

| 字段 | 含义 |
|---|---|
| `fresh_shell_mass_per_boundary_msun` | 各边界、各壳真正的新吸积质量 $A_i^{(k)}$ |
| `planned_transport_mass_per_boundary_msun` | 确定性环境质量输运矩阵 $T_{ij}^{(k)}$ |
| `transported_tracers_per_boundary` | 实际从 $i$ 迁往 $j$ 的 tracer 数量 |
| `transported_alive_weight_per_boundary` | 实际从 $i$ 迁往 $j$ 的存活人口权重 $W_{ij}^{(k)}$ |
| `transport_model` | 当前闭合名称 `nearest_shell_minimum_distance_v1` |

原有账本仍保存：

- `drawn_per_boundary`；
- `supplied_binary_weight_per_boundary`；
- `outside_merged_weight_per_boundary`；
- `entered_alive_weight_per_boundary`；
- `surviving_weight_at_boundary`；
- `raw_mergers_per_step`；
- `merged_weight_per_step`；
- `soft_gw_merged_weight_per_step`；
- `final_populations`。

注意：物理并合事件数是相应 `merged_mask` 上的权重和；`raw_mergers_per_step` 只是被抽样 tracer 的原始事件条数。

## 10. 检查点与可重复性

每完成一个完整时间边界就把以下内容写入 Drive 检查点：

- 全部人口和输运账本；
- 每壳存活者的 $(a,e,w,entry\_boundary)$；
- 下一边界编号；
- NumPy RNG 状态；
- 时间网格；
- 质量、动力学配置、抽样器参数、随机种子；
- 输运模型名、实验签名及相关源码哈希。

恢复运行时，如果配置、抽样分支、输运闭合、随机种子或源码身份变化，则拒绝读取旧检查点，避免把不同实验拼接成一个结果。

## 11. 两个 notebook 的共同设置与差异

共同设置：

- $M_0=10^{15}M_\odot$，10 壳；
- `prada12_hmf`；
- `velocity_radius_strategy="midpoint"`；
- $\xi=1$、$\mu=0.5$；
- 只运行 `full`，不再运行 `gw_only`；
- seeds 为 `20260923` 和 `20260924`；
- `N1=24576`、`N2=960`；
- 全局环境步 200 Myr；
- `adaptive_log_a`；
- 等质量 Sesana `K=sesana_endpoint`；
- 尚未考虑 soft 双星的 binary-single 环境演化。

差异只有初始轨道抽样：

| notebook | 抽样方法 |
|---|---|
| `M_010_15j_then_a.ipynb` | Fig. 19 反推锚点；先 $j$，后条件 $a$ |
| `M_010_15a_then_j.ipynb` | 严格 Appendix D 联合分布；先边缘逆 CDF 抽 $a$，后抽条件 $j$ |

两个分支必须分别解释，不能把 Fig. 19 经验复现等同为严格联合分布抽样。

## 12. 已完成的最小贯通验证

两个抽样分支都以 `N1=N2=1` 做过独立的临时贯通验证：

- 都从 $z=12$ 推进了全部 68 个演化区间；
- 都识别出 10 壳和 39 个负质量增量单元；
- 都越过了原主体函数对负增量的拒绝；
- 都生成了检查点、历史数组、速率文件、图和汇总；
- 没有正速率箱时输出明确的空图诊断，不伪造正的小数值。

这只证明控制流、账本和接口能够贯通，不是物理收敛或论文复现证明。正式高统计量结果仍需在 Colab 运行后根据云端 `summary_all.json`、`history.npz` 和速率图判断。

## 13. 当前近似、误差来源与适用条件

### 13.1 确定性模型近似

1. 壳编号距离 $|i-j|$ 只是离散层级距离，不是物理半径距离或相空间作用量距离。
2. receiver 的选择是最近壳优先工程规则，没有来自论文的校准。
3. 新吸积质量分配由“正增量减去接收的旧质量”定义，不是流体方程的解。
4. 所有同一 donor 壳内的存活 tracer 使用相同迁移概率，没有轨道—半径相关性。

### 13.2 Monte Carlo 误差

1. 实际迁移权重对确定性质量流存在有限样本波动；
2. 新 cohort 的轨道抽样存在有限样本误差；
3. 只有两个 seed，不能可靠估计完整的 seed 方差；
4. 稀有并合箱可能仍由少数高权重 tracer 主导。

### 13.3 尚未包含的物理

1. soft 双星的 binary-single 电离、破坏、弱散射和交换；
2. 连续径向轨道、壳间扩散、动力摩擦和反冲；
3. 双星再形成、再俘获、并合后代重新入晕；
4. 完整 HMF 积分和人口重加权；
5. 对输运闭合的论文或数值模拟标定。

### 13.4 当前适用条件

当前实现要求：

- halo 总壳质量在相邻边界间不下降；
- 全部 donor 需求能被正增量 receiver 容量吸收；
- 壳数在同一 $M_0$ 轨迹上固定；
- 使用 `adaptive_log_a` cohort 推进；
- 存活人口对象至少保存 $(a,e,w,entry\_boundary)$。

如果 halo 总质量下降，当前“只有外部吸积、没有整体剥离”的闭合不适用，程序应报错而不是静默处理。

## 14. 后续进入主体程序前的验收建议

在把该方案提升为 `binary_single_montecarlo.py` 的正式接口前，至少应完成：

1. 对 `N1`、`N2` 做成倍增加的收敛比较；
2. 比较不同 transport RNG seed 对分壳率和总率的影响；
3. 比较最近壳优先、累积质量坐标重映射、物理半径重叠三种闭合；
4. 检查每个边界的质量重建残差和人口守恒相对残差；
5. 检查迁移前后 $(a,e,w,cohort)$ 是否无丢失、无复制；
6. 分别检查 Fig. 19 分支和严格联合分布分支；
7. 对全局环境步做敏感性测试，确认负增量与输运量不只是 200 Myr 离散造成的伪影；
8. 明确决定未来接口是随机重映射、确定性权重拆分，还是带径向状态的动力学输运。

在这些检查完成前，本方案应继续标注为 notebook 级诊断分支，不能直接作为论文方法或最终物理模型引用。
