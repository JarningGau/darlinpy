#!/usr/bin/env python3
"""
DARLIN Python - CARLIN sequence analysis toolkit
"""

import sys

from setuptools import setup, find_packages

ext_modules = []
cmdclass = {}

try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext

    ext_modules.append(
        Pybind11Extension(
            "darlinpy.alignment._cas9_align",
            ["darlinpy/alignment/_cas9_align.cpp"],
            cxx_std=17,
        )
    )
    cmdclass["build_ext"] = build_ext
except ImportError:
    print(
        "pybind11 is not installed; building without C++ cas9_align extension.",
        file=sys.stderr,
    )

def get_version():
    with open("darlinpy/__init__.py", "r") as f:
        for line in f:
            if line.startswith("__version__"):
                return line.split("=")[1].strip().strip('"').strip("'")
    return None

def get_long_description():
    try:
        with open("README.md", "r", encoding="utf-8") as f:
            return f.read()
    except FileNotFoundError:
        return "DARLIN Python - Python implementation of CARLIN sequence analysis tools"

setup(
    name="darlinpy", 
    version=get_version(),
    author="DARLIN-toolkits Team",
    author_email="",
    description="Python implementation of CARLIN sequence analysis tools",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    url="https://github.com/jarninggau/darlinpy",
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
        "pandas>=1.3.0",
        "fuzzysearch>=0.7.0",
        "tqdm>=4.60.0",
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
            "pybind11>=2.10",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    # Ensure CARLIN configuration JSON files are installed with the package
    # The patterns are relative to the "darlinpy.config" package directory.
    package_data={
        "darlinpy.config": ["data/*.json"],
    },
    ext_modules=ext_modules,
    cmdclass=cmdclass,
)