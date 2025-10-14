#!/usr/bin/env python3
"""
DARLIN Python - CARLIN序列分析工具包
"""

from setuptools import setup, find_packages

# 读取版本号
def get_version():
    with open("darlin/__init__.py", "r") as f:
        for line in f:
            if line.startswith("__version__"):
                return line.split("=")[1].strip().strip('"').strip("'")
    return "0.1.0"

# 读取长描述
def get_long_description():
    try:
        with open("README.md", "r", encoding="utf-8") as f:
            return f.read()
    except FileNotFoundError:
        return "DARLIN Python - CARLIN序列分析工具的Python实现"

setup(
    name="darlinpy", 
    version=get_version(),
    author="DARLIN-toolkits Team",
    author_email="",
    description="CARLIN序列分析工具的Python实现",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/your-org/darlinpy",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research", 
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0", 
        "biopython>=1.79",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "isort>=5.0",
            "mypy>=0.900",
        ],
        "viz": [
            "matplotlib>=3.5.0",
            "seaborn>=0.11.0",
        ],
        "fast": [
            "numba>=0.56.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "darlin=darlin.cli:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    package_data={"darlin.config.data": ["*.json"]},
) 