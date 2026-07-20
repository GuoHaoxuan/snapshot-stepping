# 饱和恢复光变的误差模型 spec（块协方差 + 不可约涨落 + MC 闭环）

状态：**评审点已拍板**（2026-07-20，两轮无先验对抗审查后定稿）。据此改重建输出与
下游装配代码。取代旧的逐粒子 `particle_weight` 方案。

> **用途分野（先读这条）**
> - **总光变误差棒 / 对外计数比值**：下面的解析式（$U+D+\sum\text{Block}$）够用，可落地。
> - **功率谱 / QPO / 计时相位**：解析块**不够**——确定性平滑重建会在功率谱里既抹掉真
>   QPO 又造出伪 QPO，频率相关的噪声基线解析式给不准。必须走 **§9 的 Monte Carlo
>   闭环**，解析块退化为它的近似/交叉检验。

---

## 0. 一句话

恢复光变相对**真值**的完整统计协方差写成
$$\mathrm{Cov}(N)=
\underbrace{U}_{\text{不可约丢失涨落（下界）}}
+\underbrace{D}_{\text{观测计数对角泊松}}
+\sum_{\text{gap }g}\mathrm{Block}_g
+\sum_{g\neq g'}\mathrm{CoMode}_{g,g'}.$$
即"填充 bin 的泊松地板 + 普通泊松对角 + 每 gap 一个低秩块 + gap 间共模耦合"。
系统偏差（退化外推、参考盒死时间压峰、$N_{\rm lost}$ 归一、粒子≠光子）**不进协方差**，
单独标记（§8）。旧的逐粒子 $\sum w^2$ 只是其中 $D$ 加 cross-ref 块对角，且挂错了盒、
且漏了 $U$、$k$/$r$ 相关、gap 间共模。

---

## 1. 目标与范围

- 对象：**全 HE 总光变**（三盒相加）在任意时间 binning、任意能段下的均值与协方差。
  per-box 在新块方案里 **correct-by-construction**（每 gap 块自包含、明确归属目标盒），
  但主保证与测试对象是总光变。
- 覆盖（进协方差）：泊松散粒噪声、**不可约丢失涨落 $U$**、跨盒恢复引入的相关、空 bin
  插值引入的相关、标定系数 $k$ 的不确定性（含跨-ref 与 **gap 间共模**）、退化 gap 的
  率估计不确定性。
- 不覆盖（单独当系统 bias 报，§8）：退化 gap 的**系统性外推偏差**、参考盒**死时间/堆积
  压峰**、$N_{\rm lost}$ 的**归一化偏差**、**粒子≠光子**的时变分数、标定模型偏差——这些是
  bias 不是方差，任何协方差表示都装不下。
- 表示限制（§10）：暗端/能段 <1 count/bin 时高斯协方差失效，改用泊松似然。

## 2. 记号、假设与"协方差相对谁"

- 细 bin：$i$，1&nbsp;ms 重建网格。粗 bin/能段通过线性聚合（§7）。
- $C_{b,i}$：盒 $b\in\{A,B,C\}$ 在细 bin $i$ 的**观测**事件数。
- $N^{\rm fill}_{b,i}$：盒 $b$ 在细 bin $i$ 填入的 filler 数（重建计数）。
- **协方差相对真值**：$\mathrm{Cov}(N)\equiv\big\langle(N^{\rm rec}-N^{\rm true})(N^{\rm rec}-N^{\rm true})^\top\big\rangle$。
  分解 $N^{\rm rec}-N^{\rm true}=\underbrace{(N^{\rm rec}-\mathbb E[N^{\rm true}\mid\text{ref,anchor}])}_{\text{建模误差 (I)(II)(III)}}
  -\underbrace{(N^{\rm true}-\mathbb E[N^{\rm true}\mid\cdot])}_{\text{丢失事例泊松涨落 }=\,U}$。
  两部分独立 → 相加。**这是关键澄清**：填充 bin 的方差不是"filler 相对自身期望的抖动"
  一项，还有一项 $U$——真实事例被 FIFO 吞掉后、其泊松涨落**永久不可恢复**。
- **假设 P**（独立泊松）：$\{C_{b,i}\}$ 相互独立，$\mathrm{Var}(C_{b,i})=C_{b,i}$。
- 恢复量都是 $\{C_{b,i}\}$、$\{k\}$、$\{r\}$ 的**线性**函数（filler 分配 round 可忽略近似下），
  故雅可比传播即精确。
- gap $g$：目标盒 $a_g$、区间 $[t_0^g,t_1^g]$、类型 ∈ {cross-ref, degenerate}。

## 3. 均值模型（恢复光变的值）

$$N_i=\sum_b C_{b,i}\,\mathbb1[b\text{ 在 }i\text{ 未饱和}]+\sum_{g:\,i\in g}\hat F_{g,i}.$$

### 3a. cross-ref filler
目标盒 $a$、参考集 $R_g$。逐 bin shape 与插值：
$$\hat F_{g,i}=\sum_{i'}G^{g}_{i,i'}\,S^{g}_{i'},\qquad
S^{g}_{i'}=\frac{\rho_g}{n_{m,i'}}\sum_{b\in R_g}k^{g}_b\,C_{b,i'},$$
- $n_{m,i'}$：bin $i'$ 有效参考盒数；$k^g_b$：整段 gap 一个值（标定，§5b）；
  $\rho_g=N_{\rm lost}/\sum_{i'}S^{g}_{i'}\approx1$（整数守恒修正）。
- $G^g$：**插值算子**。有值 bin $G^g_{i,i}=1$；空 bin $i$ 从最近左右端点 $l,r$
  线性插值，$G^g_{i,l}=1-t,\;G^g_{i,r}=t$（$t=(i-l)/(r-l)$），单侧外推 $G^g_{i,l}=1$。

对 $\{C_{b,i'}\}$ 线性。**总灵敏度**
$$W^{g,b}_{i,i'}=\underbrace{\delta_{i,i'}\mathbb1[b\text{ meas at }i]}_{\text{自己那份}}
+\underbrace{G^{g}_{i,i'}\,\frac{\rho_g\,k^{g}_b}{n_{m,i'}}}_{\text{被 gap }g\text{ 借出的份}}.$$

### 3b. degenerate filler
无参考（三盒共饱和）或覆盖 <30%。用目标盒自己 pre/post 包率线性插值：
$$\hat F_{g,i}=\rho_g\big[r_{\rm pre}(1-t_i)+r_{\rm post}\,t_i\big]\Delta_{\rm sb},
\qquad t_i=\frac{i-i_0+\tfrac12}{n_{\rm sbins}},$$
其中 $r=n_{\rm events}/\text{span}$（**评审②：分子用该 packet 实际事例数，不用名义 109**），
$\Delta_{\rm sb}$=细 bin 宽。单侧只有一个 $r$（shape 全段常数）；(None,None) 见 §6。

## 4. 完整协方差（总公式）

四个来源：观测计数 $C$、丢失事例的固有泊松 $U$、标定 $k$、退化率 $r$。前三者 + $U$ 互相
独立，雅可比传播：
$$\boxed{\;\mathrm{Cov}(N_i,N_j)=
\underbrace{\delta_{ij}\!\sum_{g:\,i\in g}\!N^{\rm fill}_{a_g,i}}_{\text{(IV) }U\text{ 不可约丢失涨落}}
+\underbrace{\sum_{b,i'}W^{b}_{i,i'}W^{b}_{j,i'}C_{b,i'}}_{\text{(I) 泊松}}
+\underbrace{\sum_{g}\sum_{b,b'\in R_g}\!J^{k,g,b}_{i}J^{k,g,b'}_{j}\,\Sigma^{k,g}_{b,b'}}_{\text{(II) 标定}}
+\underbrace{\sum_{g\in\rm deg}\!\sum_{s}\!J^{r,g,s}_{i}J^{r,g,s}_{j}\,\mathrm{Var}(r^g_s)}_{\text{(III) 退化率}}\;}$$
其中 $W^b_{i,i'}=\delta_{i,i'}\mathbb1[b\text{ meas}]+\sum_{g:\,a_g\ne b}W^{g,b}_{i,i'}$。

**四项的角色：**
- **(IV) $U$**（新增，评审外对抗审查漏1）：每个填充 bin 贡献一个**等于其重建计数**的对角
  方差，独立于 (I)(II)(III)。物理：FIFO 吞掉的真实事例的泊松涨落不可恢复 → 方差下界，
  保证**填充 bin 误差棒 $\ge\sqrt{\text{重建计数}}$**，防止"平滑重建假装很确定→误差棒缩太小"。
  与"恢复不增信息"是一枚硬币两面：filler 的**测量噪声相关**（不该按独立 $\sqrt N$ 加，
  这是 (I) 的相关块）；但重建**相对真值**有一段独立泊松地板 $U$（我们平滑掉了、要补回）。
  $U$ 随 rebin 按泊松正常传播，但你不能靠"重建得平滑"把它压到 $\sqrt N$ 以下。
- **(I) 泊松**：观测/参考计数 $C$ 经总灵敏度 $W$ 传播。无 gap 处退化为 $\delta_{ij}C_{b,i}$
  （普通泊松 = $D$）；gap 处的插值项局域成块（§5a）。
- **(II) 标定**：$\Sigma^{k,g}_{b,b'}$ 是 gap $g$ 各参考盒 $k$ 的**满协方差**（评审①，见 §5b），
  含共用分子造成的跨-ref 相关。
- **(III) 退化率**：$\mathrm{Var}(r)=r^2/(n_{\rm events}-1)$（评审②：自由度用实际 $n_{\rm events}-1$）。

**块间共模（评审外对抗审查漏2）：** (II) 的 $k$、以及 §8 的死时间修正，在**多个 gap 间
共享**（同一对盒的相对能响、同一死时间模型）。严格地 $\mathrm{Cov}(k^g_b,k^{g'}_b)\ne0$，
产生 $\sum_{g\neq g'}\mathrm{CoMode}_{g,g'}$ 项。**对单点误差棒可忽略**（不同 gap 的 bin 不在
同一根误差棒里）；**对跨 gap 求和（fluence、总计数比值）、跨多 gap 的周期信号必须含**——
那里误差是相干相加而非块内独立。落地：先给对角块（块间独立）；跨 gap 积分时额外加
$k$ 的共模项（把同一对盒的 $k$ 当一个共享随机量，其相对不确定 $\sim1/\sqrt{C^{\rm cal}}$）。

## 5. cross-ref 块

### 5a. 泊松部分 (I)
**无插值时对角**：$G^g=I\Rightarrow W^{g,b}_{i,i'}=\delta_{i,i'}w_{b,i}$，
$$\mathrm{Cov}^{\rm P}(N_i,N_i)=\sum_b w_{b,i}^2\,C_{b,i},\quad
w_{b,i}=1+\sum_{g:\,i}\frac{\rho_g k^g_b}{n_{m,i}}.$$
这正是旧的逐粒子 $\sum w^2$——是**块的对角特例**（但旧方案漏了 $U$、(II)、(III)）。

**有插值时**：空 bin $i$ 与端点 $l,r$ 共享源计数 → gap 内非对角
$\mathrm{Cov}^{\rm P}(N_i,N_j)=\sum_b\sum_{i'}W^{g,b}_{i,i'}W^{g,b}_{j,i'}C_{b,i'}$。块大小 = gap 的 bin 数。
⚠ 插值是**确定性平滑**：它压制 gap 内高频、把功率搬向低频（§9 会咬功率谱），并把
空 bin 方差人为置零——$U$ 项把这份地板补回来。

### 5b. 标定部分 (II)，含跨-ref 满协方差（评审①：含）
$k^g_b=C_a^{\rm cal}/C_b^{\rm cal}$ 由 gap 前后 ±0.5&nbsp;s 窗现算。单个方差
$$\Big(\frac{\sigma_{k^g_b}}{k^g_b}\Big)^2=\frac1{C_a^{\rm cal}}+\frac1{C_b^{\rm cal}}.$$
**跨-ref 相关**（评审①拍板：**含**，闭式便宜）：同目标盒各 $k_b$ 共用分子 $C_a^{\rm cal}$，故
$$\Sigma^{k,g}_{b,b'}=k^g_b k^g_{b'}\Big(\frac{\delta_{bb'}}{C_b^{\rm cal}}+\frac1{C_a^{\rm cal}}\Big).$$
对角是上式，非对角 $=k_bk_{b'}/C_a^{\rm cal}$（共用分子）。雅可比
$J^{k,g,b}_i=\partial N_i/\partial k^g_b=\sum_{i'}G^g_{i,i'}\dfrac{\rho_g}{n_{m,i'}}C_{b,i'}$。
(II) 块 = $J^\top\Sigma^{k,g}J$，秩 $\le M$（$M$=参考盒数），整段 gap 全相关。实务上
$1/C_a^{\rm cal}$ 是小量（标定窗亮、计数上千），但纳入不吃亏、完整。
> **前提假设（须核实）**：$k$ 用 gap 外 ±0.5s 窗算、shape 用 gap 内参考计数，二者**时间
> 不重叠 → 独立**，故 (I) 与 (II) 无交叉项。实作里确认标定窗与 shape 参考 bin 不相交。

## 6. 退化块 (III)（评审②③拍板）

雅可比 $J^{r,g,\rm pre}_i=\rho_g(1-t_i)\Delta_{\rm sb}$，$J^{r,g,\rm post}_i=\rho_g t_i\Delta_{\rm sb}$。
$\mathrm{Var}(r)=r^2/(n_{\rm events}-1)$。整段 gap 由 $(r_{\rm pre},r_{\rm post})$ 驱动 → **秩-2** 块。
- **评审②（拍板）**：$r=n_{\rm events}/\text{span}$ 的分子与自由度都用该 packet **实际事例数**
  $n_{\rm events}$，不用名义 109（reset 前后常是残包）。改 `degenerate_gap_variance` 与
  `reconstruct_gaps` 里的 `EVENTS_PER_PKT` 硬编码。
- 单侧（只一个 $r$）：秩-1，$J_i=\rho_g\Delta_{\rm sb}$（shape 常数）。
- **(None,None)（评审③拍板）**：无相邻包率，率**无从估**。→ **不入协方差**（不编造方差），
  在 §8 给该段 `sys_bias_scale`≈100% 的系统标记，并给 filler 打**可屏蔽**标志，下游做
  计时/功率谱可直接 mask。

## 7. 下游装配（任意 bin/能段）

粗 bin $I$（可含能段 $e$）：均值 $N_I=\sum_{i\in I}N_i$，协方差
$\mathrm{Cov}(N_I,N_J)=\sum_{i\in I}\sum_{j\in J}\mathrm{Cov}(N_i,N_j)$。实作：
1. **均值** $N_I$ = 在 $I$、能段 $e$ 内**数所有事件**（EVT + FILL_GAP）。
2. **泊松地板**：观测事件（EVT）计数 → $D$；填充事件（FILL_GAP）计数 → $U$。二者都是
   "数该 bin 该能段的行数"，对角。
3. **各 gap 块**：对与 $I$（及 $J$）相交的 gap，用块描述子 + 事件流重算源计数 $C_{b,i}$，
   按 §5/§6 组 (I 插值项)+(II)+(III)，累加进 $\mathrm{Cov}$。
4. **跨 gap 积分/比值**：额外加 §4 的 $k$ 共模项。
5. **误差棒**：$\sigma_I=\sqrt{\mathrm{Cov}(N_I,N_I)}$。**功率谱/χ²**：见 §9/§10。

### 7a. 能段扩展（评审⑤拍板：eband 确定性，无抽样噪声）
`assign_gap_fill_channels` 是**确定性分位重采样**（取参考 in-gap 分布等间隔分位、van der
Corput 位反转铺时间，**无 RNG**）。故：
- **总光变**：filler 不管能量都被数 → 能量赋值**零方差贡献**。上面全套原样成立。
- **能段光变**（只在做能段/硬度产物时才建，当前非主保证）：
  1. **方差放大**：分位映射把 $N_{\rm ref}$ 个参考铺成 $N_{\rm filler}$ 个 filler，一个参考事例被
     复制到多个 filler 时其泊松涨落放大 $\sim N_{\rm filler}/N_{\rm ref}$ 倍进该能段。块里的
     源计数 $C_{b,i}$ 换成 $C_{b,i,e}$，并带这个复制放大因子。
  2. **能段间反相关**：总 filler 数固定、按参考谱确定性拆分 → 能段间负相关。硬度比/谱
     延迟必须建**跨能段协方差块**（$\mathrm{Cov}(N_{I,e},N_{I,e'})<0$）。
  3. **系统偏**：filler 谱 = 参考 in-gap 谱（或退化时侧窗谱），若目标盒 gap 内真实谱不同
     → 能段 filler 系统偏，进 §8。

## 8. 系统偏差（不入协方差，单独报告）

天体物理惯例：stat + sys **分列不相加**。每 gap 输出 `sys_bias_flag`（bool）+
`sys_bias_scale`（量级估计）。**评审④拍板**：`sys_bias_scale` 现在给一个便宜的粗代理
$\approx|r_{\rm post}-r_{\rm pre}|/(r_{\rm pre}+r_{\rm post})$（率变化越大→越弯→线性外推误差越大，
**下界启发式**，非严格上界）；**精确量化**用 §11 的"掩掉-重建注入"实验，留待科学分析逐 case。

系统项清单（对抗审查加重）：
- **退化线性外推**（可主导退化 gap）：burst 是曲线、填直线 → 整段系统偏（凸峰必被削）。
  三盒共饱和只在**最亮、曲率最大**处启用，削峰最重。
- **参考盒死时间/堆积压峰**（对抗审查漏，cross-ref 主系统项之一）：未 FIFO 复位 ≠ 线性
  响应；峰上"未饱和"参考盒本身死时间/pileup 压平 → 拿它当模板系统性低估目标盒峰。
  gap 前 ±0.5s 的 $k$ 窗、退化的 $r_{\rm pre/post}$ 锚点同样已在非线性区（"gap 外=干净"是错觉）。
- **$N_{\rm lost}$ 归一化偏差**（乘性，直接进对外比值）：$N_{\rm lost}$ 是**推断值**不是观测值；
  饱和期包率是 MCU 节流地板 → $N_{\rm lost}$ 系统偏小 → 恢复计数偏低。它错，整条填充光变
  按比例错。是否含本底粒子、filler 计不计入比值分母（f7 已剔除 filler）都要一致声明。
- **粒子≠光子的时变分数**：重建的是**粒子**率；转光子需一个随时间变、峰上未知的粒子/光子
  比。声明产物是"计数级"而非"光子级"，或给该系统带。
- **$k$ 常数假设 & pre/post 跨台阶**：$k$ 用 gap 前后窗算并外推进 gap，源谱在暴期演化 →
  $k(t)$ 变；gap 恰在最亮最快变处，pre/post 两侧不同流强平均成一个 $k$ 会拉平峰。
- **能标 / 增益下沉 / NaI-CsI 甄别**（仅能段）：高率 gain sag 挪能段边界；脉冲甄别误分。
- **比值估计 Jensen 偏**：$E[C_a/C_b]$ 偏高 $\sim1/C_b^{\rm cal}$（小）。
- **(None,None) 地板猜测**：纯猜，≈100%（评审③）。
- **off-axis 角响应**：HE 窄视场准直，GRB 常大离轴，逐盒角响应差异使盒间比非常数
  （记忆里"从 CsI 一侧到达"已否证，须显式当系统）。

## 9. 功率谱 / QPO 用途：Monte Carlo 闭环（解析块不够）

确定性平滑重建在功率谱里**既能抹掉真 QPO 又能造出伪 QPO**，且频率相关噪声基线解析式
给不准。凡涉及功率谱/QPO/周期性，**以 MC 闭环为主结果，解析块仅作交叉检验**：

- **频率相关噪声基线 $P_{\rm noise}(f)$**：纯泊松无信号数据过**完整重建管线** → 得经验白噪声
  基线。**放弃 Leahy 平坦 $P=2$ 假设**（这直接回应内部"Leahy 高频白噪声待查证"：混合段
  基线依赖频率与 gap 占比，插值平滑把高频压到 2 以下）。一切显著性对 $P_{\rm noise}(f)$ 判定。
- **伪周期自查**：packet 准周期到达 → packet 频率及谐波的梳状伪 QPO；1ms binning 的 sinc
  边瓣；多 reset → 多 gap 转移函数在 gap 间距差频处拍出假峰。每个宣称 QPO 须给"**离源/空区
  同法重建**"的对照谱，伪峰复现即证伪；并列 packet 到达间隔分布证明与宣称频率不重合。
- **检测效率（转移函数）**：注入-回收——多频率注入已知振幅 QPO，给检测效率 vs 频率；任何
  "未检测到 QPO"的上限须除以该频率效率（平滑压真信号 → 检测效率非平、有假阴性）。
- **有效独立频点数 / 归一**：填充 bin 不独立 → 独立频点 $<N/2$；Leahy 归一的光子数**只用
  真实探测事件**，不含 filler。误差棒用 MC 协方差，不用 $\chi^2_2$。
- **计时相位轴**（协方差之外的独立一根轴）：filler 事件时间是**内插放的**，不是真到达时刻
  → 重建区脉冲相位/峰到达时刻有系统抖动。**这条直接关系 221009A/250919A 的时间解算验证**：
  用验证暴给出重建区计时偏差的实测上界，作为独立不确定度轴单独报，不塞进计数协方差。

## 10. 表示方法限制：低计数用泊松似然

1ms bin 即使 $10^4$ cts/s 也才 ~10 counts/bin；能段/暗端 <1 count/bin，分布强偏、离散。
高斯协方差（哪怕带完整非对角）在低计数区表示不了偏态离散分布。规则：
- **亮段/粗 bin（每 bin ≫ 10 counts）**：高斯 + §4 协方差，可用误差棒/χ²。
- **暗端/能段/细 bin（每 bin ≲ few counts）**：用**泊松似然**（Cash/C-stat），不用高斯误差棒；
  或证明高斯近似在该产物上误差 <X%。

## 11. 验证实验清单（证明模型数值上没漏）

误差**枚举**再全也要**闭环验证**。落地顺序：
1. **null 闭环**：纯泊松无信号 → 过管线 → 时域比真值、频域给 $P_{\rm noise}(f)$。
2. **掩掉-重建注入**：取一个**未饱和**已知亮暴，人为掩掉目标盒，重建后逐 bin 比真值
   （时域 + 功率谱）。**唯一能直接量化 $U$（重建相对真值误差）和退化外推 bias 的实验。**
3. **经验协方差 vs 模型**：MC 多次实现算经验 $C_{ij}$，与解析块比特征谱，匹配才算数。
4.（做功率谱时）QPO 注入-回收给检测效率；5. 有效自由度 / 独立频点数。

## 12. 方法参数敏感性（研究者自由度）

以下旋钮的结论敏感性须表征，否则算未量化系统误差：
- 标定窗宽 ±0.5s、cross-ref↔degenerate 的 30% 覆盖阈：各 ±50% 扫描，展示总光变/比值/
  （将来）QPO 结论稳定；阈值附近 gap 归类跳变 → 方法依赖不连续。
- **telemetry 空档 vs 饱和 gap**：SAA 关机/模式切换/遥测丢包的空档**不得**当饱和 gap 去填。
  重建前须确认 gap 确由 FIFO reset 造成。

## 13. 输出格式

**(a) 事件流**（均值）——去掉 `particle_weight` 列：
```
box, type, met, channel, pulse_width, pkt_idx, evt_idx
```
（`type` ∈ {EVT, FILL_GAP}；均值=数行，$U$/$D$ 泊松地板都由数行得到。）

**(b) gap 块表**（协方差）——新文件，每 gap 一行：
```
gap_id, target_box, type,            # crossref | degenerate
t_start, t_stop,
ref_boxes,                           # crossref: 参与参考盒
k[ref], C_a_cal, C_ref_cal[ref],     # crossref: k 与算 Σ_k 所需分子/分母计数
r_pre, r_post, n_pre, n_post,        # degenerate: 端点率与实际事例数(评审②)
maskable,                            # (None,None) 等不可信段的可屏蔽标志(评审③)
sys_bias_flag, sys_bias_scale        # 系统项标记与量级(评审④)
```
- 源计数 $C_{b,i}$ 不存，下游把事件流在 $[t_0,t_1]$ 内按盒分箱重算。
- 插值算子 $G^g$ 由"哪些 bin 空"（无有效参考事件）重建，不必显存。
- $\Sigma^{k,g}$ 由 $k$、$C_a^{\rm cal}$、$C_{\rm ref}^{\rm cal}$ 现算（评审①）。

**(c) 下游库**：函数 (事件流, 块表, 目标 binning/能段) → (均值向量, 协方差稀疏表示：
$U$/$D$ 对角 + 每 gap 块 + 可选 gap 间 $k$ 共模)。功率谱产物改走 §9 的 MC 管线。

## 14. 与旧方案的关系

- 旧 `particle_weight`：EVT 存 $w$、filler 存 $v$，称"方差=$\sum$weight²"。它等于本 spec 的
  **(I) 泊松对角**，但漏了 **(IV) $U$**、(I) 插值非对角、(II) 标定（及跨-ref/gap 间共模）、
  (III) 退化 gap 内相关；且把 cross-ref 方差挂在**参考盒**事例上，单盒会错。
- 本 spec：观测事件回归权重 1（普通泊松 $D$），恢复引入的 (协)方差收进每 gap 自包含块，
  补上填充 bin 的泊松地板 $U$ → 单盒正确、相关完整、更紧凑（块=每 gap 几个数）。功率谱
  用途明确 pivot 到 MC 闭环。

---

## 附：评审点拍板记录（2026-07-20）

1. **§5b 跨-ref 相关** → **含**。$\Sigma^{k,g}_{b,b'}=k_bk_{b'}(\delta_{bb'}/C_b^{\rm cal}+1/C_a^{\rm cal})$，
   闭式便宜，(II) 升为 $J^\top\Sigma^k J$。实务小但完整。
2. **§6 fencepost** → 一律用**实际 $n_{\rm events}$**（不用名义 109）；$\mathrm{Var}(r)=r^2/(n_{\rm events}-1)$。
3. **§6 (None,None)** → **不入协方差** + 系统标记（≈100%）+ 可屏蔽标志。
4. **§8 `sys_bias_scale`** → 现给 bool + 粗代理 $|r_{\rm post}-r_{\rm pre}|/(r_{\rm pre}+r_{\rm post})$；
   精确量化用 §11 掩掉-重建注入，留待科学分析。
5. **§7a 能段** → eband 确定性、**无抽样噪声**；总光变零影响；能段产物另建（方差放大 +
   能段间反相关 + 谱系统偏）三项。
6. **范围** → 主保证 = 全 HE 总光变（已定协方差）；per-box correct-by-construction，测试待补。

## 附：对抗审查补进的结构性项（两轮无先验子 agent，2026-07-20）

- **漏1 → (IV) $U$ 不可约丢失涨落**：填充 bin 泊松地板，独立、下界，防误差棒缩太小。
- **漏2 → gap 间共模**：$k$/死时间跨 gap 共享，跨 gap 积分/周期信号须含（§4 CoMode）。
- **漏3 → 低计数泊松似然**（§10）：暗端高斯协方差失效。
- **漏4 → 功率谱 MC 闭环**（§9）：$P_{\rm noise}(f)$、伪周期、检测效率、独立频点、计时相位轴。
