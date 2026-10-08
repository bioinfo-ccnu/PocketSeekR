# RNA-only pocket identification: figure legend

## English caption

**Figure X | Workflow of RNA-only pocket identification.** **a**, The input comprises the observed RNA heavy-atom coordinates and atom identities, without a ligand. **b**, A 1.5 Å grid samples probe-accessible free space around the RNA. Directional ray tests characterize local enclosure, and connected region growth generates geometrical candidates. Geometric scoring and redundancy filtering retain up to 40 regions for chemical evaluation. **c**, Aromatic-face and directional polar support are evaluated across the candidate region. Local co-support combines nearby complementary fields through spatial neighborhoods connected by probe-clear segments. The regional quality score is \(Q(r)=G(r)[1+2c(r)]\), where \(G(r)\) is the geometric score and \(c(r)\) is the bounded regional chemical score. **d**, Greedy joint selection balances quality and novelty of the RNA context crops using \(\operatorname{gain}(r\mid S)=Q(r)[1-J_{\max}(r,S)]^{0.25}\). Here, \(J_{\max}\) is the maximum Jaccard overlap between the candidate crop and previously selected crops; additional distance and overlap constraints exclude redundant candidates. **e**, The method returns up to five ranked pocket centers with RNA context crops formed by a 10 Å neighborhood and completion of the observed nucleotides.

Molecular structures, grid points and support distributions are schematic illustrations rather than coordinates or results from a particular experimental structure. Translucent halos indicate RNA context neighborhoods rather than cavity boundaries. Chemical fields represent potential interaction support and do not depict a bound ligand or establish an interaction energy. Parameters were selected through offline calibration; inference uses fixed heuristic rules without a trained predictive model.

## 中文图注

**图 X｜仅基于 RNA 的口袋识别流程。** **a**，输入为 RNA 已观测重原子的坐标和原子身份，不输入配体。**b**，利用 1.5 Å 网格采样 RNA 周围探针可达的自由空间，通过方向射线衡量局部包围程度，再通过连通区域生长形成候选。几何评分和冗余过滤后，保留最多 40 个区域进行化学评价。**c**，在整个候选区域上计算芳香碱基面支持和方向性极性支持，并通过探针可通行的局部邻域计算协同支持。区域质量评分为 \(Q(r)=G(r)[1+2c(r)]\)，其中 \(G(r)\) 为几何评分，\(c(r)\) 为有界的区域化学评分。**d**，联合选择使用贪心增益 \(\operatorname{gain}(r\mid S)=Q(r)[1-J_{\max}(r,S)]^{0.25}\)，兼顾质量和 RNA 裁剪区域的新颖程度；其中 \(J_{\max}\) 为候选裁剪与已选裁剪之间的最大 Jaccard 重叠，另有距离与重叠约束去除冗余候选。**e**，输出最多 5 个排序后的口袋中心，以及由 10 Å 邻域并补全已观测核苷酸构成的 RNA 上下文裁剪。

图中的分子、网格和化学支持分布均为示意，不对应某个真实结构的实验结果。半透明范围表示 RNA 上下文邻域，不代表空腔的边界。化学支持表示潜在相互作用条件，不表示已有结合配体或实际相互作用能。参数经过离线校准；推断过程使用固定启发式规则，不依赖训练得到的预测模型。

## Files and generation

- Image: `workflow.png`
- Original raster size: 1672 × 941 pixels.
- Prompt: `prompt.txt`
- Generation mode: built-in `image_gen.imagegen` tool, with the `imagegen` skill; no CLI generation.
- Algorithm source: `docs/pocket_v5_algorithms.tex`.
- Image and captions were checked for the five-stage flow, the 40-region/5-output upper limits, the two score formulas and absence of a ligand input.
