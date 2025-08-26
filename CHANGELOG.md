# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2024-01-XX

### Added
- Initial release of DARLIN Python
- Core CARLIN sequence alignment functionality
- Position-specific gap penalty system
- NUC44 and multiple scoring matrices support
- CARLIN amplicon configuration system
- Sequence sanitization and standardization
- AlignedSEQ and AlignedSEQMotif data structures
- Comprehensive test suite
- Batch sequence processing
- Detailed alignment statistics and analysis

### Core Features
- **cas9_align algorithm**: High-precision sequence alignment with position-specific penalties
- **AmpliconConfig**: Complete CARLIN amplicon configuration management
- **ScoringConfig**: Multiple scoring matrices (NUC44, simple, transition/transversion)
- **CARLINAligner**: Integrated high-level alignment interface
- **SequenceSanitizer**: Two-stage sequence standardization
  - Prefix/postfix contamination removal
  - Conserved region error correction
- **AlignedSEQ system**: Motif-level sequence analysis and event classification

### Technical Implementation
- Position-specific gap penalties based on biological motifs
- Automatic motif boundary calculation
- Robust error handling and input validation
- Comprehensive documentation and examples
- Development environment setup scripts

### Testing and Quality
- Unit tests for all core components
- Integration tests for end-to-end workflows
- Sanitization functionality tests
- Performance benchmarks
- Code examples and demonstrations

### Documentation
- Complete API documentation
- Usage examples and tutorials
- Development setup guide
- Biological background and methodology

### Compatibility
- Python 3.8+
- NumPy 1.20.0+
- SciPy 1.7.0+
- BioPython 1.79+

### Performance
- Efficient NumPy-based algorithms
- Optimized motif boundary calculations
- Batch processing capabilities
- Memory-efficient sequence handling 