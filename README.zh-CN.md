# PocketSeekR

**从 RNA 三维结构识别小分子结合口袋。**

[![CI](https://github.com/bioinfo-ccnu/PocketSeekR/actions/workflows/ci.yml/badge.svg)](https://github.com/bioinfo-ccnu/PocketSeekR/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/bioinfo-ccnu/PocketSeekR?color=416788)](https://github.com/bioinfo-ccnu/PocketSeekR/releases/latest)
[![Python](https://img.shields.io/badge/Python-%3E%3D3.10-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![pip](https://img.shields.io/badge/pip-installable-3776AB?logo=pypi&logoColor=white)](#安装)
[![CPU](https://img.shields.io/badge/CPU-supported-52796F)](#运行)
[![Input](https://img.shields.io/badge/Input-RNA%20structure-52796F)](#输入与输出)

[English](README.md) | 简体中文

PocketSeekR 是独立的 Python 口袋识别工具。它从 RNA 的探针可达空间生成候选区域，计算整个区域内的芳香、方向性极性及局部交叉支持，再结合评分与区域重叠，从最多 **40 个内部候选**中选出至多 **5 个口袋**。预测使用 CPU，不需要输入配体或训练模型。

## 安装

通过 pip 安装 GitHub 上的固定版本：

```bash
python -m pip install "git+https://github.com/bioinfo-ccnu/PocketSeekR.git@v0.1.0"
pocketseekr --version
```

也可从 [Releases](https://github.com/bioinfo-ccnu/PocketSeekR/releases/latest) 下载 wheel 文件，使用 `python -m pip install pocketseekr-0.1.0-py3-none-any.whl` 安装。要求 Python ≥3.10，依赖 NumPy、SciPy、PyYAML 和 PyTorch；PyTorch 仅用于坐标张量，不需要 CUDA 或 RDKit。Linux 集群可先安装 CPU 版 PyTorch：

```bash
python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
```

## 运行

```bash
pocketseekr --rna-pdb path/to/rna.pdb --output-dir outputs/my_rna
```

也可使用 `python -m pocketseekr` 调用相同入口。克隆仓库后可直接运行随附的真实 RNA 示例：

```bash
git clone https://github.com/bioinfo-ccnu/PocketSeekR.git
cd PocketSeekR
pocketseekr --rna-pdb examples/1F1T_rna.pdb --output-dir outputs/1F1T
```

## Python 调用

```python
from pocketseekr import predict

result = predict("path/to/rna.pdb")
for pocket in result["pockets"]:
    print(pocket["rank"], pocket["center"], pocket["score"])
```

默认参数随安装包提供，不依赖当前工作目录。需要显式指定参数时，可用 `--config path/to/settings.yaml`；字段格式见 [默认配置](configs/default.yaml)。Python 接口也支持 `from pocketseekr import PocketSettings`，再将 `PocketSettings()` 传入 `predict(..., settings=...)`。

## 输入与输出

输入为 PDB 格式的 RNA 结构，采用首个模型。支持的核苷酸名称为 `A`、`U`、`G`、`C`、`I`、`RA`、`RU`、`RG`、`RC`、`PSU` 和 `H2U`；只使用这些核苷酸的观测重原子。输入中的配体及其他分子不参与口袋检测。mmCIF 结构需先转换为合法的 PDB 文件。

| 输出文件 | 内容 |
| --- | --- |
| `pockets.json` | 最终口袋的中心、排名、评分、RNA 裁剪原子身份、参数及来源哈希。 |
| `candidates.json` | 内部区域候选池与最终选择结果。 |

中心坐标单位为 Å，沿用输入坐标系；RNA 上下文按中心周围 10 Å 裁剪，并扩展到完整观测核苷酸。JSON 原子索引对应解析后的 RNA 原子顺序，裁剪集合不等同于空腔边界。

## 方法与开发

方法包括探针可达区域构建、整个区域的化学支持评价，以及固定预算下的联合选择，详见 [方法说明](docs/method.md)。几何与化学质量按 `Q(r) = G(r) × [1 + 2c(r)]` 组合，再以裁剪重叠抑制重复输出。

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
python -m build
```

源码从 RiboPoseDiff 的当前口袋预测路径独立提取，来源与数值一致性验证见 [来源说明](docs/provenance.md)。仓库只提供 Python 包、命令行和 API，不包含网页功能、对接训练或历史评估数据。
