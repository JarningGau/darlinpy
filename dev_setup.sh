#!/bin/bash
# DARLIN Python 开发环境快速启动脚本

echo "🚀 设置DARLIN Python开发环境..."

# 检查Python版本
python_version=$(python3 --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1,2)
echo "📋 Python版本: $python_version"

# 创建虚拟环境 (如果不存在)
if [ ! -d "venv" ]; then
    echo "📦 创建虚拟环境..."
    python3 -m venv venv
fi

# 激活虚拟环境
echo "🔄 激活虚拟环境..."
source venv/bin/activate

# 升级pip
echo "⬆️  升级pip..."
pip install --upgrade pip

# 安装依赖
echo "📚 安装项目依赖..."
pip install -r requirements.txt

# 以开发模式安装包
echo "🔧 安装DARLIN Python包 (开发模式)..."
pip install -e .

# 运行测试验证
echo "🧪 运行基础测试..."
python tests/test_cas9_align.py

echo ""
echo "✅ 开发环境设置完成！"
echo ""
echo "🎯 使用方法:"
echo "   source venv/bin/activate  # 激活环境"
echo "   python tests/test_cas9_align.py  # 运行测试"
echo "   python -c \"import darlin; print('DARLIN Python已准备就绪!')\"  # 验证安装"
echo ""
echo "📁 项目结构:"
echo "   darlin/                   # 主包"
echo "   ├── alignment/            # 序列比对模块"
echo "   │   └── cas9_align.py    # ✅ 核心比对算法"
echo "   ├── config/              # 配置模块"
echo "   ├── calling/             # 等位基因调用模块"
echo "   ├── mutations/           # 突变注释模块"
echo "   └── utils/               # 工具模块"
echo ""
echo "🎉 准备开始开发吧！" 