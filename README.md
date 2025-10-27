# DARLIN Python

This repository is based on https://gitlab.com/hormozlab/carlin and https://github.com/ShouWenWang-Lab/Custom_CARLIN. The original CARLIN pipeline is implemented in MATLAB. Here, we provide a Python implementation of CARLIN sequence analysis tools for CRISPR-Cas9 lineage tracing. This repository is only focus on DARLIN sequence alignment and mutation calling.


## Features

- 🧬 **CRISPR-Cas9 Sequence Alignment**: High-precision alignment algorithm with position-specific gap penalties
- 🎯 **Allele Calling**: Robust allele identification based on event patterns
- 📝 **Mutation Annotation**: HGVS-format mutation event annotation
- 🔧 **Multiple Arrays Support**: Built-in support for Col1a1, Rosa, and Tigre arrays. And can extend to self-defined array.

## Quick Start

### Installation

```bash
# Install from GitHub (recommended)
git clone https://github.com/jarninggau/darlinpy.git
cd darlinpy
pip install -e .

# Or install with development dependencies
pip install -e ".[dev]"

# Verify installation
python -c "import darlin; print(f'DARLIN Python v{darlin.__version__} installed successfully!')"
```

### Basic Usage

```python
from darlin import analyze_sequences

# Analyze CARLIN sequences
sequences = [
    "CGCCGGACTGCACGACAGTCGACCGATGGAGTCGACACGACTCGCGCATATTCGATGGAGTCGACTACAGTCGCTACGAGTATGGAGTCGATACGTAGCACGCAGAACGATGGGAGCT",
    "CGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGATACGCGCACGCTATGGAGTCGAGAGCGCGCTCGTCGACTATGGAGTCGCGACTGTACGCACACGCGATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGAGTCGAGACGCTGACGATATGGAGTCGATACGTAGCACGCAGACGATGGGAGCT"
]

# One-line analysis with default Col1a1 configuration
results = analyze_sequences(sequences, config='Col1a1', verbose=True)

print(f"Analysis completed: {len(results.valid_sequences)} sequences processed")
print(f"Average alignment score: {results.alignment_scores.mean():.2f}")
print(f"Total mutations detected: {len(results.mutations)}")
```

### Advanced Usage

```python
from darlin import analyze_sequences, AmpliconConfig, build_carlin_config

# Custom amplicon configuration
## Col1a1 array (CA)
## segments = conserved site (13bp) + cut site (7bp)
## pam = PAM (TGG) + linker (4bp)
## the reference is prefix + (segments+pam)x10 + postfix
## the sequencing library contains:
## primer5 + prefix + (segments+pam)x10 + postfix + secondary sequence + primer3
sequence = {
    "segments" : [
        "GACTGCACGACAGTCGACGA",
        "GACACGACTCGCGCATACGA",
        "GACTACAGTCGCTACGACGA",
        "GCGAGCGCTATGAGCGACTA",
        "GATACGATACGCGCACGCTA",
        "GAGAGCGCGCTCGTCGACTA",
        "GCGACTGTACGCACACGCGA",
        "GATAGTATGCGTACACGCGA",
        "GAGTCGAGACGCTGACGATA",
        "GATACGTAGCACGCAGACGA"
        ],
    "pam" : "TGGAGTC",
    "prefix" : "CGCCG",
    "postfix" : "TGGGAGCT",
    "Primer5" : "GAGCTGTACAAGTAAGCGGC",
    "Primer3" : "CGACTGTGCCTTCTAGTTGC",
    "SecondarySequence" : "AGAATTCTAACTAGAGCTCGCTGATCAGCCT"
}
build_carlin_config(sequence, output_path="configs/array_Col1a1.json")

# Build and use custom configuration
build_carlin_config(sequence_config, output_path="custom_array.json")
config = AmpliconConfig(config_file="custom_array.json")

# Analyze with custom configuration
results = analyze_sequences(sequences, config=config, method='exact')
```

## Supported Lineage Arrays

DARLIN Python supports multiple array configurations:

- **Col1a1 (CA)**: 276 bp reference with 10 segments
- **Rosa (RA)**: 275 bp reference with 10 segments  
- **Tigre (TA)**: 275 bp reference with 10 segments

## License

MIT License

## Citation

If you use this tool, please cite the paper:

1. L. Li#, S. Bowling#, H. Lin, D. Chen, S.-W. Wang*, F. D. Camargo*, DARLIN mouse for in vivo lineage tracing at high efficiency and clonal diversity. Nature Protocols, doi: 10.1038/s41596-025-01141-z (2025).

2. L. Li, S. Bowling, S. E. McGeary, Q. Yu, B. Lemke, K. Alcedo, Y. Jia, X. Liu, M. Ferreira, A. M. Klein, S.-W. Wang*, F. D. Camargo*, A mouse model with high clonal barcode diversity for joint lineage, transcriptomic, and epigenomic profiling in single cells. Cell 186, 5183-5199.e22 (2023). 

3. S. Bowling*, D. Sritharan*, F. G. Osorio, M. Nguyen, P. Cheung, 
A. Rodriguez-Fraticelli, S. Patel, W-C. Yuan, Y. Fujiwara, B. E. Li, S. H. Orkin, 
S. Hormoz, F. D. Camargo. "An Engineered CRISPR-Cas9 Mouse Line for 
Simultaneous Readout of Lineage Histories and Gene Expression Profiles 
in Single Cells." Cell (2020), https://doi.org/10.1016/j.cell.2020.04.048 