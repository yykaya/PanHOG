# PanHOG Configuration Guide

## Overview

PanHOG now supports YAML configuration files to simplify parameter management and make it easier to run complex analyses with consistent settings.

## Key Features

- **YAML Configuration**: Store all your analysis parameters in a structured YAML file
- **Command-line Override**: Command-line arguments take precedence over config file settings
- **Flexible Usage**: Use config files for defaults while still allowing per-run customization
- **Error Handling**: Graceful handling of missing or malformed config files

## Usage

### Basic Usage with Config File

```bash
# Use default config.yaml file
python PanHOG.py --hog data/N0.tsv --fasta data/fasta_files/

# Use custom config file
python PanHOG.py --config my_custom_config.yaml --hog data/N0.tsv --fasta data/fasta_files/
```

### Command-line Override

Even when using a config file, you can override any setting via command-line:

```bash
# Override bootstrap count and output directory
python PanHOG.py --config config.yaml --hog data/N0.tsv --fasta data/fasta_files/ --bootstrap 200 --output ./new_results/
```

## Configuration File Structure

The config file uses YAML format. Here's an example with all available options:

```yaml
# Input files (still required on command line or here)
hog: "path/to/your/N0.tsv"
fasta: "path/to/your/fasta_directory"

# Output settings
output: "./results"
prefix: "my_analysis_"

# Analysis options
pan: true
clade: "species1,species2,species3"  # Optional: for clade-specific analysis

# Optional analyses
proteome: "ALL"  # or specify species list
genevar: "ALL"   # or specify species list
saturation: true
saturation_cladepair: false

# Bootstrap settings
bootstrap: 10000

# Visualization settings
marker_core: "o"
marker_pan: "s"

# Clade-specific settings
clade1: "species1,species2,species3"
clade2: "species4,species5,species6"

# Colors and markers for cladepair saturation analysis
marker_core_clade1: "^"
marker_pan_clade1: "^"
marker_core_clade2: "s"
marker_pan_clade2: "s"
color_core_clade1: "#c0392b"
color_pan_clade1: "#f1c40f"
color_core_clade2: "#c0392b"
color_pan_clade2: "#3498db"

# Transformation options
zscore: false
log: true
```

## Configuration Priority

Settings are applied in this order (later overrides earlier):

1. **Default values** (built into the script)
2. **Config file values** (from YAML)
3. **Command-line arguments** (highest priority)

## Example Workflows

### 1. Standard Analysis with Config

Create a `config.yaml`:
```yaml
output: "./my_results"
prefix: "experiment1_"
saturation: true
genevar: "ALL"
bootstrap: 2000
```

Run analysis:
```bash
python PanHOG.py --hog data/N0.tsv --fasta data/fasta_files/
```

### 2. Multiple Experiments with Different Configs

Create different config files for different experiments:

`experiment1.yaml`:
```yaml
output: "./exp1_results"
prefix: "exp1_"
bootstrap: 14000
saturation: true
```

`experiment2.yaml`:
```yaml
output: "./exp2_results"
prefix: "exp2_"
bootstrap: 50000
saturation_cladepair: true
```

Run experiments:
```bash
python PanHOG.py --config experiment1.yaml --hog data/N0.tsv --fasta data/fasta_files/
python PanHOG.py --config experiment2.yaml --hog data/N0.tsv --fasta data/fasta_files/
```

### 3. Config with Command-line Overrides

Use a base config but override specific parameters:
```bash
# Use config defaults but change output directory and bootstrap count
python PanHOG.py --config base_config.yaml --hog data/N0.tsv --fasta data/fasta_files/ --output ./special_run --bootstrap 10000
```

## Error Handling

- **Missing config file**: Script continues with command-line arguments only (warning displayed)
- **Malformed YAML**: Script exits with error message
- **Invalid parameters**: Standard argparse validation applies

## Benefits

1. **Reproducibility**: Save exact analysis parameters for future reference
2. **Batch Processing**: Easy to run multiple analyses with different settings
3. **Complex Configurations**: Manage many parameters without long command lines
4. **Documentation**: Config files serve as documentation of analysis settings
5. **Flexibility**: Still allows command-line overrides when needed

## Migration from Original Script

Your existing command-line workflows will continue to work unchanged. The config file feature is completely optional and backward-compatible.

## Dependencies

The modified script requires the `pyyaml` package:

```bash
pip install pyyaml
```

All other dependencies remain the same as the original script.
