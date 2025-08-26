# DARLIN Python

CARLIN序列分析工具的Python实现

## 功能特点

- 🧬 **CRISPR-Cas9序列比对**: 位置特异性gap惩罚的高精度比对算法
- 🎯 **等位基因调用**: 基于事件模式的鲁棒等位基因识别
- 📝 **突变注释**: HGVS格式的突变事件注释
- ⚡ **高性能**: 基于NumPy优化的算法实现

## 快速开始

### 安装

```bash
# 从GitHub安装 (推荐)
git clone https://github.com/your-org/darlinpy.git
cd darlinpy
pip install -e .

# 或安装开发依赖
pip install -e ".[dev]"

# 验证安装
python -c "import darlin; print(f'DARLIN Python v{darlin.__version__} 安装成功!')"
```

### 基本使用

```python
from darlin.alignment import align_to_carlin

# 比对单个CARLIN序列
sequence = "CGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATAC..."
result = align_to_carlin(sequence, verbose=True)

print(f"比对得分: {result['alignment_score']:.2f}")
print(f"序列一致性: {result['statistics']['identity']*100:.1f}%")
```

### 批量序列比对

```python
from darlin.alignment import CARLINAligner

# 创建比对器
aligner = CARLINAligner()

# 批量比对
sequences = [
    "CGCCGGACTGCACGACAGTCGACGATGGAGTC...",
    "CGCCGGACTGCACGACAGTCGACGATGGAGTC...",
]

results = aligner.align_sequences(sequences)
for i, result in enumerate(results):
    print(f"序列 {i+1}: 得分 {result['alignment_score']:.1f}")
```

### 底层API使用

```python
from darlin.alignment.cas9_align import cas9_align, nt2int, int2nt

# 直接使用比对算法
seq = nt2int("ACGTACGT")
ref = nt2int("ACGTGCGT")

score, al_seq, al_ref = cas9_align(seq, ref, open_penalty, close_penalty, sub_score)
print(f"比对得分: {score:.2f}")
print(f"比对结果: {int2nt(al_seq)} vs {int2nt(al_ref)}")
```

## 项目状态

🔄 **开发中** - 核心比对算法已完成，正在实现完整的分析流程

### 已完成
- ✅ 核心cas9_align算法
- ✅ 项目结构搭建
- ✅ 基础配置文件
- ✅ CARLIN扩增子配置系统
- ✅ NUC44评分矩阵
- ✅ 位置特异性gap惩罚
- ✅ 集成的CARLIN比对器
- ✅ 批量序列处理
- ✅ 详细比对统计

### 开发中
- 🔄 序列标准化
- 🔄 等位基因调用
- 🔄 突变注释

## 开发

```bash
# 克隆仓库
git clone <repository>
cd darlinpy

# 安装开发环境
pip install -e ".[dev]"

# 运行测试
pytest

# 代码格式化
black darlin/
isort darlin/
```

## 许可证

MIT License

## 引用

如果使用本工具，请引用原始CARLIN论文：

> S. Bowling*, D. Sritharan*, F. G. Osorio, M. Nguyen, P. Cheung, 
A. Rodriguez-Fraticelli, S. Patel, W-C. Yuan, Y. Fujiwara, B. E. Li, S. H. Orkin, 
S. Hormoz, F. D. Camargo. "An Engineered CRISPR-Cas9 Mouse Line for 
Simultaneous Readout of Lineage Histories and Gene Expression Profiles 
in Single Cells." Cell (2020), https://doi.org/10.1016/j.cell.2020.04.048 