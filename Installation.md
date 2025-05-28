# PanHOG - Phylogeny-informed Pangenome Analysis Tool

Welcome to PanHOG, a powerful tool for analyzing pan-genomes and their hierarchical orthologous groups (HOGs).

## Documentation

- [Configuration Guide](README_config.md)
- [Pangene Integration Guide](README_pangene.md)

## Setup and Usage

1. First, create the conda environment:
```bash
conda env create -f environment_pangene.yml
conda activate pangene
```

2. For detailed usage instructions, please refer to the documentation links above.

## Project Structure

- `PanHOG.py`: Main PanHOG analysis tool
- `PangeneHOG.py`: Integration module for pangene analysis
- `environment_pangene.yml`: Conda environment configuration
- `example*config.yaml`: Configuration templates
