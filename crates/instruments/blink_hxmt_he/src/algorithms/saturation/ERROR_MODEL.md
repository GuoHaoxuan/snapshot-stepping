# 饱和恢复光变的误差模型 spec（三张表 + 1ms 网格解析 S·diag(C)·Sᵀ）

状态：**评审点已拍板**（2026-07-20 两轮无先验对抗审查定稿；2026-07-21 重构为
"三张表 + 1ms 网格解析协方差"）。取代旧的 `D + Σ块 + U 对角` 方案，也取代更旧的
逐粒子 `particle_weight`。据此改重建输出（三张表）与下游装配代码（1ms 网格拼 S）。

> **一句话方针（先读这条）**
> 恢复光变 $N$ 是独立源 $x=(C,k,r,U)$ 的**线性像**。一切协方差都从 **1ms 完整网格**
> 上的稀疏灵敏度矩阵 $S$ 解析算出：$\mathrm{Cov}(N)=S\,\mathrm{diag}(C)\,S^\top+(\text{k 项})
> +(\text{r 项})+(\text{U 对角，可开关})$。任意 binning（含 1ms、0.1s、1s、乃至傅里叶域功率谱）
> 一律 $P\,\mathrm{Cov}\,P^\top$，$P$ 是把 1ms 格聚合到目标基的线性算子。
> **纯解析、确定、可复现——全程禁用蒙特卡洛 / RNG**（§10）。
> 系统偏差（退化外推、死时间压峰、$N_{\rm lost}$ 归一、粒子≠光子）与协方差**无任何数学
> 通路**，走严格平行台账单列（§9）。

---

## 0. 为什么是 S·diag(C)·Sᵀ（核心洞见）

三盒 FIFO 饱和后，target 盒 gap 段的计数是**借**参考盒计数、经标定 $k$ 与插值 $G$
线性重建出来的。旧方案（对角 $D$ + 每 gap 低秩块 + 对角 $U$）的病：把 filler 的方差
挂到**被借的参考盒事例**上，而不是 filler 实际**落点**（target 盒、gap 时间）；且用块
对角近似丢了 filler↔参考的**完全相关**。后果：三盒相加时正负相关没对齐，总光变方差偏小
约 30%，注入验证 pull≈1.32（误差棒系统性偏窄）。

修法：把整个重建写成一个**稀疏线性映射** $N=S\,x$。方差记在 $S$ 的**行**（=落点）、
源记在 $S$ 的**列**。filler↔参考的完全相关是 $S$ 里的**跨盒非对角**，自动出现、符号正确。
三盒相加 = 对行做盒求和，每 bin 方差退化为**逐粒子 $\sum w^2 C$**（数值正确，pull→1）。

---

## 1. 目标与范围

- 对象：**全 HE 总光变**（三盒相加）在任意时间 binning、任意能段下的均值与协方差；
  per-box 在 $S$ 表示下 **correct-by-construction**（每行显式归属落点盒）。
- 覆盖（进协方差 $S\,\mathrm{diag}(C)\,S^\top+k+r+U$）：泊松散粒噪声、跨盒恢复相关、
  空 bin 插值相关、标定 $k$ 不确定（含跨-ref 与 gap 间共模）、退化 gap 率估计不确定、
  不可约丢失涨落 $U$（可开关）。
- 不覆盖（单独当系统 bias 报，§9）：退化线性外推偏差、参考盒死时间/堆积压峰、
  $N_{\rm lost}$ 归一化偏差、粒子≠光子时变分数、$k$ 常数假设跨台阶、能标下沉、off-axis
  角响应——这些是 **bias 不是 variance**，与协方差无数学通路。
- 表示限制（§11）：暗端/能段 <1 count/bin 时高斯协方差失效，改泊松似然。

## 2. 记号、1ms 网格、"协方差相对谁"

- **1ms 完整网格**：细 bin $i$，$\Delta_{\rm sb}=1$&nbsp;ms（`SHAPE_BIN_WIDTH`）。一切先在此网格
  精确表示，粗 bin/能段/功率谱都由聚合算子 $P$（§8）从这里线性得到。
- $C_{b,i}$：盒 $b\in\{A,B,C\}$ 在 1ms bin $i$ 的**观测**事件数（从事件流数），泊松 $\mathrm{Var}=C_{b,i}$。
- $N^{\rm fill}_{b,i}$：盒 $b$ 在 bin $i$ 填入的 filler 数（重建计数）。
- **协方差相对真值**：$\mathrm{Cov}(N)\equiv\langle(N^{\rm rec}-N^{\rm true})(N^{\rm rec}-N^{\rm true})^\top\rangle$。
  分解 $N^{\rm rec}-N^{\rm true}=(N^{\rm rec}-\mathbb E[N^{\rm true}\mid\text{ref,anchor}])-\underbrace{(N^{\rm true}-\mathbb E[\cdot])}_{=\,U}$。
  前者是 $S\,x$ 的建模误差（$C,k,r$ 传播），后者是被 FIFO 吞掉事例的**永久不可恢复**泊松
  涨落 $U$。二者独立 → 相加。$U$ 做成**可开关独立腿**（§6d）。
- **假设 P（独立泊松）**：$\{C_{b,i}\}$ 相互独立，$\mathrm{Var}(C_{b,i})=C_{b,i}$。
- **线性**：恢复量都是 $x=(\{C_{b,i}\},\{k\},\{r\},\{U\})$ 的线性函数（filler 分配 round 可忽略
  近似下），雅可比传播即精确。
- gap $g$：目标盒 $a_g$、区间 $[t_0^g,t_1^g]$、类型 ∈ {crossref, degenerate}。

---

## 3. 三张表 schema（权威：Rust 产出、Python 消费必须逐字一致）

均值由**事件流**得到；协方差由**gap 参数表** + **gap 格结构表**在 1ms 网格上拼 $S$ 得到。
`reconstruct_gaps` 内部已算好 $n_m$、哪些格空、插值端点/τ（`bin_refs` 与
`interpolate_empty_bins`），现在把它们**序列化出来**，别让 Python 重推（重推=与 Rust 漂移，
同"轨道文件硬依赖"同类坑）。

### 表① 事件流 `events.csv`（现有，不变）
```
box,type,met,channel,pulse_width,pkt_idx,evt_idx
```
- `type` ∈ {EVT, FILL_GAP}。
- 给：均值（数行）；源计数 $C_{b,i'}$（按盒×binning 分箱数 EVT 行）；$S$ 的**恒等对角**
  （每个 EVT 在其 $(b,i)$ 贡献 $S[(b,i),(b,i)]=1$）。

### 表② gap 参数表 `gapcov.csv`（现有基础上**加 `rho` 列**）
每 gap 一行，变长字段（`ref_boxes`/`k`/`c_ref_cal`）用分号分隔：
```
gap_id,target_box,type,t_start,t_stop,ref_boxes,k,c_ref_cal,c_a_cal,rho,
r_pre,r_post,n_pre,n_post,maskable,sys_bias_flag,sys_bias_scale
```
- `type` ∈ {crossref, degenerate}。
- crossref 段：填 `ref_boxes`（参与形状重建的参考盒名，分号列表）、`k`（各盒标定系数，
  分号列表）、`c_ref_cal`（各盒标定窗计数，分号列表）、`c_a_cal`（target 标定窗计数，
  所有 ref 共用分子）、`rho`（$\rho_g=N_{\rm lost}/\sum S$ 整数守恒修正，$\approx1$）；退化字段留空。
- degenerate 段：填 `r_pre,r_post`（端点率）、`n_pre,n_post`（端点 packet 实际事例数=自由度）、
  `maskable`（(None,None) 地板段可屏蔽）；crossref 字段留空。
- `sys_bias_flag`=（degenerate ? true : false）；`sys_bias_scale` 见 §9（评审④粗代理）。
- 给：k 项、r 项、U、系统台账的参数。

### 表③ gap 格结构表 `gapbins.csv`（**新增**）
每 (gap, 1ms 格) 一行，**每格都记**（不做"只记非默认"优化，下游据此精确拼 $S$，不猜插值）：
```
gap_id,bin_index,t_lo,n_m,kind,left_bin,right_bin,tau
```
- `kind` ∈ {measured, empty}。
- **measured**（该格有有效参考）：`n_m`=该格有效参考盒数（=`bin_refs[si].len()`）；
  `left_bin,right_bin,tau` 留空。
- **empty**（无参考，插值填）：`left_bin,right_bin`=最近左/右端点格的 `bin_index`；
  `tau`=插值权重 $\tau=(i-l)/(r-l)$；`n_m` 留空（下游用端点格的 `n_m`）。单侧外推 → 只有
  一端非空、$\tau$ 记 0 或 1、缺的一端 `bin_index` 留空。
- 给：$S$ 的 filler↔参考**逐格系数** $G^g/n_m$（下游在 1ms 网格精确拼 $S$）。

> **一致性铁律**：$n_m$、`left/right/tau` 是 Rust 已经算过的量（`bin_refs` 给 $n_m$ 与哪些格
> 空，`interpolate_empty_bins` 给端点/τ），序列化即可。Python **绝不**重算插值算子。

---

## 4. 均值模型（恢复光变的值）

$$N_i=\sum_b C_{b,i}\,\mathbb1[b\text{ 在 }i\text{ 未饱和}]+\sum_{g:\,a_g,\,i\in g}\hat F_{g,i}.$$

### 4a. crossref filler
目标盒 $a$、参考集 $R_g$。逐 bin 形状 + 插值：
$$\hat F_{g,i}=\sum_{i'}G^{g}_{i,i'}\,S^{g}_{i'},\qquad
S^{g}_{i'}=\frac{\rho_g}{n_{m,i'}}\sum_{b\in R_g}k^{g}_b\,C_{b,i'}.$$
- $n_{m,i'}$：bin $i'$ 有效参考盒数（表③ `n_m`）；$k^g_b$：整段 gap 一个值（表② `k`）；
  $\rho_g$：整数守恒修正（表② `rho`）。
- $G^g$：**插值算子**。measured 格 $G^g_{i,i}=1$；empty 格 $i$ 从端点 $l,r$ 线性插值
  $G^g_{i,l}=1-\tau,\;G^g_{i,r}=\tau$（表③ `tau`），单侧外推 $G^g_{i,l}=1$。

### 4b. degenerate filler
无参考或覆盖 <30%。用目标盒自己 pre/post 包率线性插值：
$$\hat F_{g,i}=\rho_g\big[r_{\rm pre}(1-t_i)+r_{\rm post}\,t_i\big]\Delta_{\rm sb},\qquad
t_i=\frac{i-i_0+\tfrac12}{n_{\rm sbins}},$$
$r=n_{\rm events}/\text{span}$（评审②：分子用该 packet **实际事例数**，不用名义 109）。单侧只有
一个 $r$（shape 全段常数）；(None,None) → 不入协方差，`maskable`（§6c/§9）。

---

## 5. 稀疏灵敏度矩阵 S 的构造（权威）

$S$ 行 = 落点 $(b,\,\text{1ms bin } i\,[,\text{能段 } e])$；列 = 源计数 $C_{b,i'}$。三类非零元：

**(i) 恒等**（每个观测 EVT）：
$$S[(b,i),(b,i)] \mathrel{+}= 1.$$

**(ii) crossref filler，measured 格**（表③ `kind=measured`，$i'=i$，$G=$ 恒等）：对 $b\in R_g$
$$S[(a,i),(b,i)] \mathrel{+}= \frac{\rho_g\,k^g_b}{n_{m,i}}.$$

**(iii) crossref filler，empty 格**（表③ `kind=empty`，从端点 $l,r$、$\tau$ 内插）：对 $b\in R_g$
$$S[(a,i),(b,l)] \mathrel{+}= (1-\tau)\frac{\rho_g\,k^g_b}{n_{m,l}},\qquad
S[(a,i),(b,r)] \mathrel{+}= \tau\,\frac{\rho_g\,k^g_b}{n_{m,r}}.$$
（$n_{m,l},n_{m,r}$ 取端点 measured 格的 `n_m`。单侧外推只留 $l$ 端、系数 $\rho_g k^g_b/n_{m,l}$。）

**要点：**
- **filler 方差落在 target 行 $(a,i)$**：$\mathrm{Var}(N_{a,i})^{\rm P}=\sum_{b,i'}S[(a,i),(b,i')]^2\,C_{b,i'}$。
  方差记在**落点**（target 盒、gap 时间），不再借挂参考盒位置——修好旧方案的归位病。
- **filler↔参考完全相关 = 跨盒非对角**：同一源计数 $C_{b,i'}$ 既进参考行 $(b,i')$（系数 1）
  又进 target 行 $(a,i)$（系数 $S[(a,i),(b,i')]$），故
  $\mathrm{Cov}(N_{a,i},N_{b,i'})^{\rm P}=S[(a,i),(b,i')]\cdot1\cdot C_{b,i'}\neq0$。**自动出现、符号正确**，
  保住完全相关，总光变不再偏小。
- **总光变（三盒相加）= 对行做盒求和**：每 bin 方差
  $$\mathrm{Var}\Big(\sum_b N_{b,i}\Big)^{\rm P}=\sum_b\Big(1+\sum_{g:\,a_g\ne b}\frac{\rho_g k^g_b}{n_{m,i}}\Big)^2 C_{b,i}
  =\sum_b w_{b,i}^2\,C_{b,i}.$$
  这**恰是逐粒子 $\sum w^2 C$**——是 $S\,\mathrm{diag}(C)\,S^\top$ 的行-盒求和特例，数值正确
  （修复旧 D+U 偏小 ~30%、pull 1.32→1 的病）。measured 格 $w_{b,i}=1+\sum_g\rho_g k^g_b/n_{m,i}$；
  empty 格 $w$ 由端点系数按 $\tau$ 混合。

---

## 6. 完整协方差（总公式）

四条**独立**腿，雅可比传播直接相加：
$$\boxed{\;\mathrm{Cov}(N)=
\underbrace{S\,\mathrm{diag}(C)\,S^\top}_{\text{(I) 观测计数泊松}}
+\underbrace{\textstyle\sum_g J^{k,g}\,\Sigma^{k,g}\,{J^{k,g}}^\top}_{\text{(II) 标定 }k}
+\underbrace{\textstyle\sum_{g\in\rm deg} J^{r,g}\,\mathrm{diag}(\mathrm{Var}\,r)\,{J^{r,g}}^\top}_{\text{(III) 退化率 }r}
+\underbrace{\lambda_U\,\mathrm{diag}(N^{\rm fill})}_{\text{(IV) }U\text{ 可开关}}\;}$$

### 6a. (I) 观测计数泊松 $S\,\mathrm{diag}(C)\,S^\top$
$C$ 从事件流（表①）在 1ms 网格数；$S$ 从表①（恒等）+ 表②（$k,\rho$）+ 表③（$n_m$/端点/τ）
按 §5 拼。无 gap 处退化为 $\delta_{ij}C_{b,i}$（普通泊松对角）；gap 内 empty 格与端点共享源计数
→ 局域非对角块。总光变行-盒求和 = §5 的 $\sum w^2 C$。

### 6b. (II) 标定 $k$，含跨-ref 满协方差（评审①：含）
$k^g_b=C_a^{\rm cal}/C_b^{\rm cal}$ 由 gap 前后 ±0.5&nbsp;s 窗现算（表② `k,c_a_cal,c_ref_cal`）。满协方差
$$\Sigma^{k,g}_{b,b'}=k^g_b k^g_{b'}\Big(\frac{\delta_{bb'}}{C_b^{\rm cal}}+\frac1{C_a^{\rm cal}}\Big),$$
对角 $=k_b^2(1/C_b^{\rm cal}+1/C_a^{\rm cal})$，非对角 $=k_bk_{b'}/C_a^{\rm cal}$（共用分子 $C_a^{\rm cal}$
造成跨-ref 相关）。雅可比 $J^{k,g,b}_i=\partial N_i/\partial k^g_b=\sum_{i'}G^g_{i,i'}\dfrac{\rho_g}{n_{m,i'}}C_{b,i'}$
（同 $S$ 里去掉 $k_b$ 的那份）。(II) 块 $=J^{k,g}\Sigma^{k,g}{J^{k,g}}^\top$，秩 $\le|R_g|$，整段 gap 全相关。
> **独立性前提（须核实）**：$k$ 用 gap 外 ±0.5s 窗、shape 用 gap 内参考计数，时间不重叠 → 独立，
> 故 (I) 与 (II) 无交叉项。实作确认标定窗与 shape 参考 bin 不相交。
>
> **gap 间共模**：同一对盒的 $k$ 在多个 gap 间共享，严格 $\mathrm{Cov}(k^g_b,k^{g'}_b)\neq0$。
> 对**单点误差棒可忽略**（不同 gap 的 bin 不在同一根误差棒里）；对**跨 gap 求和**（fluence、
> 总计数比值）、周期信号必须含——那里 $k$ 相干相加。落地：先给块内独立，跨 gap 积分时把同一
> 对盒的 $k$ 当一个共享随机量（相对不确定 $\sim1/\sqrt{C^{\rm cal}}$）额外加共模项。

### 6c. (III) 退化率 $r$
秩-2 块（评审②③）：每 gap
$$\sigma^2_{\rm gap}=\Big(\frac T2\Big)^2\Big(\frac{r_{\rm pre}^2}{n_{\rm pre}-1}+\frac{r_{\rm post}^2}{n_{\rm post}-1}\Big),$$
雅可比 $J^{r,g,\rm pre}_i=\rho_g(1-t_i)\Delta_{\rm sb}$，$J^{r,g,\rm post}_i=\rho_g t_i\Delta_{\rm sb}$，
$\mathrm{Var}(r_s)=r_s^2/(n_s-1)$。方差**分到 gap 内 filler 行**（局域，不畸变相邻真实事例）。
- 评审②：分子/自由度都用**实际 $n_{\rm events}$**（不用名义 109）。
- 单侧（只一个 $r$）：秩-1，$J_i=\rho_g\Delta_{\rm sb}$。
- **(None,None)（评审③）**：率无从估 → **不入协方差**（不编造方差），`maskable=true` + §9 系统
  标记 $\approx100\%$；下游计时/功率谱可直接 mask。

### 6d. (IV) $U$ 不可约丢失涨落（**可开关独立腿**）
每个 filler 格贡献一个**等于其重建计数**的对角方差 $\mathrm{diag}(N^{\rm fill})$（泊松地板）。
物理：FIFO 吞掉的真实事例的泊松涨落**永久不可恢复**。开关 $\lambda_U\in\{0,1\}$：
- **$\lambda_U=0$（默认，测量方差）**：画误差棒 / χ² / 对外计数比值。此时协方差是"重建量相对
  其自身期望"的**测量 scatter**，filler 的测量噪声已在 (I) 的相关块里（不重复计）。
- **$\lambda_U=1$（恢复 vs 真值）**：注入验证 / 回应审稿人"重建有多准"。补回被平滑掉的独立泊松
  地板，保证 filler 格误差棒 $\ge\sqrt{\text{重建计数}}$，防"平滑重建假装很确定"。

输出务必标注：$U$ 是**非测量 scatter**，是重建-真值差的下界，不是数据点的测量误差。

---

## 7. 方差归位（行=落点，列=源）

一句话总纲，防回退旧病：
- **行 = 落点**：方差记在 filler 实际填补的位置 $(a_g,\,\text{gap 时间}\,i)$，**不借挂**参考盒位置。
- **列 = 源**：不确定性来源是参考盒的观测计数 $C_{b,i'}$（列），一个源经 $S$ 同时进它自己的
  参考行和所有借它的 target 行。
- **跨盒非对角 = 完全相关**：$S$ 的非零 off-diagonal 自动、正确地保住 filler↔参考的完全相关；
  三盒相加时正负相关精确对齐，总光变不再偏小。
- 旧 `particle_weight` 把 cross-ref 方差挂在**参考盒事例**上、块对角近似丢相关、漏 $U$/$k$/$r$
  —— per-box 错、总光变偏小。本 spec 用 $S$ 归位一次性修好。

---

## 8. 任意 binning：聚合算子 P（含功率谱）

粗 bin / 能段 / 任意基一律
$$\mathrm{Cov}(N^{P})=P\,\mathrm{Cov}(N)\,P^\top,$$
$P$ 把 1ms 格线性聚合到目标基。**1ms / 0.1s / 1s 全从同一 1ms 完整表示算**，不各自近似。
- 时间 rebin：$P$ = 0/1 求和矩阵（把细格并进粗格），能段选择同理（按能段过滤行）。
- **功率谱 / 傅里叶域**：DFT 也是线性算子——取 $P=F$（DFT 矩阵），则频域噪声协方差
  $F\,\mathrm{Cov}(N)\,F^\top$ **也是纯解析**，无需 MC。频率相关噪声基线、有效独立频点数
  由此直接读出（filler 不独立 → 独立频点 $<N/2$；Leahy 归一光子数只用真实 EVT，不含 filler）。
- 唯一进不了 $P\,\mathrm{Cov}\,P^\top$ 的是**系统偏差**（平滑重建抹真信号/造伪结构、计时相位抖动）
  —— 那是 bias 不是 variance，走 §9 平行台账。

## 9. 系统偏差：严格平行台账（无数学通路）

天体物理惯例 stat + sys **分列不相加**。系统偏差与协方差**无任何数学通路**，单列。每 gap 输出
`sys_bias_flag`（bool）+ `sys_bias_scale`（量级）。**评审④**：`sys_bias_scale` 现给粗代理
$\approx|r_{\rm post}-r_{\rm pre}|/(r_{\rm pre}+r_{\rm post})$（下界启发式，非严格上界）；精确量化用 §12
掩掉-重建注入，留待科学分析逐 case。

系统项清单：
- **退化线性外推**（可主导退化 gap）：burst 是曲线、填直线 → 凸峰必被削；三盒共饱和只在最亮、
  曲率最大处启用，削峰最重。
- **参考盒死时间/堆积压峰**（cross-ref 主系统项）：峰上"未饱和"参考盒本身死时间/pileup 压平 →
  当模板系统性低估 target 峰；$k$ 窗、$r_{\rm pre/post}$ 锚点同样已在非线性区（"gap 外=干净"是错觉）。
- **$N_{\rm lost}$ 归一化偏差**（乘性，直接进对外比值）：$N_{\rm lost}$ 是推断值；饱和期包率是 MCU
  节流地板 → $N_{\rm lost}$ 系统偏小 → 恢复计数偏低。filler 计不计入比值分母须一致声明（f7 已剔 filler）。
- **粒子≠光子时变分数**：重建的是**粒子**率；转光子需峰上未知的时变粒子/光子比。声明产物是
  "计数级"而非"光子级"。术语用 "particle" 不用 "photon"。
- **$k$ 常数假设 & pre/post 跨台阶**：源谱暴期演化 → $k(t)$ 变；gap 恰在最亮最快变处，两侧平均
  成一个 $k$ 会拉平峰。
- **平滑重建的功率谱偏 & 计时相位轴**：确定性插值压 gap 内高频、把功率搬向低频（抹真 QPO / 造伪
  结构）；filler 事件时间是内插放的、非真到达时刻 → 重建区脉冲相位/峰到达时刻有系统抖动。**这条
  直接关系 250919A 的时间解算验证**：用验证暴给出重建区计时偏差的实测上界，作为独立不确定度轴
  单独报，**不塞进计数协方差**。
- **能标 / 增益下沉 / NaI-CsI 甄别**（仅能段）；**off-axis 角响应**（盒间比非常数，"从 CsI 一侧
  到达"已否证，显式当系统）；**(None,None) 地板猜测** ≈100%（评审③）。

## 10. 纯解析、确定、可复现——**删除蒙特卡洛**

**全程禁用 MC / RNG。** 恢复是 $x\mapsto S\,x$ 的线性映射，其协方差有闭式 $S\,\mathrm{diag}(C)\,S^\top+k+r+U$；
任意 binning（含功率谱，§8 的 $P=F$）都是同一闭式的线性像。因此：
- **不做**任何"多次实现取经验协方差"的 MC 闭环——它引入 RNG、不可复现、且是解析式已给的量的
  有噪估计。
- 旧 spec 的 §9 "功率谱 MC 闭环" **作废**：频率相关噪声基线由 $F\,\mathrm{Cov}\,F^\top$ 解析给出，
  伪周期/检测效率由 $S$ 的转移函数解析分析（$S$ 的行谱即转移函数），无需注入采样。
- eband 赋值是**确定性分位重采样**（等间隔分位 + van der Corput 位反转铺时间，**无 RNG**，评审⑤）；
  总光变里 filler 不管能量都被数 → 能量赋值零方差贡献。能段产物（另建）带方差放大 + 能段间反相关
  + 谱系统偏三项，仍全解析。
- 验证注入（§12）也是**确定性**的：固定掩码、固定重建，逐 bin 比真值，无随机采样。

可复现判据：同一输入（事件流 + 轨道文件）→ 逐位相同的三张表 → 逐位相同的协方差。

## 11. 表示限制：低计数用泊松似然

1ms bin 即使 $10^4$ cts/s 也才 ~10 counts/bin；能段/暗端 <1 count/bin，分布强偏、离散。
高斯协方差（哪怕带完整非对角）表示不了偏态离散分布。规则：
- **亮段/粗 bin（≫10 counts）**：高斯 + §6 协方差，误差棒/χ² 可用。
- **暗端/能段/细 bin（≲few counts）**：用**泊松似然**（Cash/C-stat），或证高斯近似误差 <X%。

## 12. 验证实验清单（确定性注入，证明数值上没漏）

误差枚举再全也要闭环验证。全部**确定性**（无 MC，§10）：
1. **掩掉-重建注入**（主）：取未饱和已知亮暴，**确定性**掩掉 target 盒，重建后逐 bin 比真值（时域）。
   已量化：$U$（$\lambda_U=1$ 腿）+ **cross-ref measured 腿** + **共饱和 empty 腿**（`cosaturate`：把
   参考盒也标 unreliable，造真机制 empty）。实测 250919A 平坦 baseline（138 gap）：measured
   pull.std=1.06、共饱和 empty 0.998、bias~0.1%，可从仓库复现（`scripts/injection_validation.sh`）。
   ⚠**未覆盖：degenerate（pre/post 外推）腿。** `cosaturate` 造的是 cross-ref empty，不是真退化
   （需 $\mathrm{has\_ref}$=false 的全宽共饱和），故退化 r 项与其外推 bias **尚未被真值验证**；
   补验前须先修 `inject` 的 `prev/next_pkt_idx=0` 硬编码（否则退化注入用文件首包率）。
   记忆铁律：**验证只用 250919A**（211211A 已弃用：峰上 59% 三盒共饱和退化、SPI-ACS 拿不到）。
2. **解析协方差自洽**：$S\,\mathrm{diag}(C)\,S^\top$ 行-盒求和应逐位等于逐粒子 $\sum w^2 C$；
   pull 分布（注入残差 / 解析 σ）均值 0、宽度 1（旧方案 1.32，本 spec 目标 →1）。
3. **1ms↔0.1s↔1s 一致性**：$P\,\mathrm{Cov}\,P^\top$ 从 1ms 聚合 vs 直接在粗网格算，须一致。
4. **null 自洽**：纯泊松无信号过管线 → 频域读 $P_{\rm noise}(f)=[F\,\mathrm{Cov}\,F^\top]_{ff}$（解析），
   与经验周期图比（此步允许经验图仅作交叉核对，非结果来源）。

## 13. 方法参数敏感性（研究者自由度）

以下旋钮结论敏感性须表征，否则算未量化系统误差：
- 标定窗宽 ±0.5s、cross-ref↔degenerate 的 30% 覆盖阈：各 ±50% 扫描，展示总光变/比值稳定；
  阈值附近 gap 归类跳变 → 方法依赖不连续。
- **telemetry 空档 vs 饱和 gap**：SAA 关机/模式切换/遥测丢包空档**不得**当饱和 gap 去填；重建前
  确认 gap 确由 FIFO reset 造成。**绝不设 `MAX_SEC_GAP`**（该变量让 blink 重建结果错误）。
- 浮点**绝不判相等**，用 $(a-b)$ 绝对值 $<$ 容差。

---

## 附：评审点拍板记录（2026-07-20）

1. **§6b 跨-ref 相关** → **含**。$\Sigma^{k,g}_{b,b'}=k_bk_{b'}(\delta_{bb'}/C_b^{\rm cal}+1/C_a^{\rm cal})$，
   闭式便宜，(II) = $J^\top\Sigma^k J$。
2. **§6c fencepost** → 一律用**实际 $n_{\rm events}$**；$\mathrm{Var}(r)=r^2/(n_{\rm events}-1)$。
3. **§6c (None,None)** → **不入协方差** + 系统标记（≈100%）+ `maskable`。
4. **§9 `sys_bias_scale`** → bool + 粗代理 $|r_{\rm post}-r_{\rm pre}|/(r_{\rm pre}+r_{\rm post})$；精确量化用
   §12 掩掉-重建注入。
5. **§10 eband** → 确定性、**无抽样噪声**；总光变零影响；能段产物另建（方差放大 + 能段间反相关 +
   谱系统偏）。
6. **范围** → 主保证 = 全 HE 总光变；per-box 在 $S$ 表示下 correct-by-construction。

## 附：对抗审查补进的结构性项（两轮无先验子 agent，2026-07-20）

- **漏1 → (IV) $U$ 不可约丢失涨落**：filler 泊松地板，独立、下界；本 spec 做成**可开关腿**（§6d）。
- **漏2 → gap 间 $k$ 共模**：跨 gap 积分/周期信号须含（§6b）。
- **漏3 → 低计数泊松似然**（§11）：暗端高斯协方差失效。
- **漏4 → 功率谱**：由 $P=F$ 的解析 $F\,\mathrm{Cov}\,F^\top$ 处理（§8），**取代**旧 MC 闭环（§10 作废之）。

## 附：重构记录（2026-07-21：D+Σ块+U → 三张表 + S·diag(C)·Sᵀ）

- **动机**：旧对角/低秩块方案把 filler 方差挂错位置（参考盒而非落点）、块对角丢 filler↔参考完全
  相关 → 总光变偏小 ~30%、pull 1.32。
- **修法**：整个重建写成稀疏线性映射 $N=S\,x$，方差归位（行=落点、列=源，§7），跨盒非对角自动
  保住完全相关，总光变行-盒求和退化为逐粒子 $\sum w^2 C$（数值正确，pull→1）。
- **落地**：三张表（§3）——事件流不变、gapcov 加 `rho`、新增 gapbins（每 1ms 格的 $n_m$/端点/τ，
  由 `bin_refs`/`interpolate_empty_bins` 序列化，Python 不重推）；下游在 1ms 网格拼 $S$，任意
  binning 走 $P\,\mathrm{Cov}\,P^\top$（§8）。
- **MC**：删除（§10），纯解析确定可复现。
