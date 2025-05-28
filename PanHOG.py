#!/usr/bin/env python3
import csv
import sys
import glob
import os
import argparse
import random
import statistics
from collections import defaultdict
import yaml

import matplotlib
matplotlib.use('Agg')  #environments without an X server
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

##############################
# Configuration Loading
##############################

def load_config(config_file):
    """
    Load configuration from YAML file.
    Returns a dictionary with configuration parameters.
    """
    if not os.path.exists(config_file):
        print(f"[WARNING] Config file '{config_file}' not found. Using command-line arguments only.")
        return {}
    
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        print(f"[INFO] Loaded configuration from '{config_file}'")
        return config if config else {}
    except yaml.YAMLError as e:
        print(f"[ERROR] Failed to parse YAML config file '{config_file}': {e}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Failed to read config file '{config_file}': {e}")
        sys.exit(1)

def merge_config_with_args(config, args):
    """
    Merge configuration from YAML with command-line arguments.
    Command-line arguments take precedence over config file.
    """
    # Only update args with config values if the arg wasn't explicitly set on command line
    for key, value in config.items():
        if hasattr(args, key):
            # Check if the argument was set to its default value (meaning it wasn't specified on command line)
            current_value = getattr(args, key)
            
            # For boolean flags, only update if current value is False (default)
            if isinstance(value, bool) and not current_value:
                setattr(args, key, value)
            # For other types, check if it's still the default value
            elif not isinstance(value, bool) and current_value in [None, False, ".", "", 100]:  # common defaults
                setattr(args, key, value)
    
    return args

##############################
# Classes for storing HOG data
##############################

class HOG:
    def __init__(self, ID):
        self.ID = ID
        self.subgenomes = []
        self.cultivars = set()
        self.members = defaultdict(dict)
        self.passport = []
        
    def addSubGenome(self, s):
        self.subgenomes.append(s)
    def addOGID(self, o):
        self.OGID = o
    def addClade(self, clade):
        self.clade = clade
    def addCultivar(self, c):
        self.cultivars.add(c)
    def addMembers(self, k, v, n):
        self.members[k][v] = n
    def addPassport(self, p):
        self.passport.append(p)

#############################################
# FASTA parsing: extract protein sequences
#############################################

def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return a dict: {gene_id: sequence}
    """
    seqs = {}
    with open(fasta_file) as f:
        current_id = None
        current_seq = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = "".join(current_seq)
    return seqs

###########################################################
# getMissingGenes: collects "cloud" genes by comparing FASTA
# to HOG membership (no name cleaning).
###########################################################

def getMissingGenes(hogsfile, fasta_dir, clade_filter=None):
    dSpecies = {}
    dSpeciesGenes = defaultdict(list)
    dAllGenes = defaultdict(list)
    dGenesMissing = defaultdict(list)

    with open(hogsfile, "r") as fdh:
        reader = csv.reader(fdh, delimiter="\t")
        header = next(reader)
        # columns 0,1,2 => HOG, OG, Clade
        for i in range(3, len(header)):
            sp_name = header[i]
            if clade_filter and sp_name not in clade_filter:
                continue
            dSpecies[i] = sp_name

        for line in reader:
            for i in range(3, len(line)):
                if i in dSpecies and line[i] != '':
                    genes = [g.strip() for g in line[i].split(", ")]
                    dSpeciesGenes[dSpecies[i]].extend(genes)

    fasta_files = glob.glob(os.path.join(fasta_dir, "*.fa")) + glob.glob(os.path.join(fasta_dir, "*.fasta"))
    for fasta in fasta_files:
        sp = os.path.basename(fasta).split('.')[0]
        if clade_filter and sp not in clade_filter:
            continue
        seq_dict = parse_fasta(fasta)
        dAllGenes[sp] = list(seq_dict.keys())

    for idx in dSpecies:
        sp_name = dSpecies[idx]
        if sp_name not in dAllGenes:
            sys.stderr.write(f"[WARNING] No FASTA found for species '{sp_name}'. Skipping missing gene calculation.\n")
            continue
        missing = set(dAllGenes[sp_name]) - set(dSpeciesGenes[sp_name])
        dGenesMissing[sp_name] = missing

    return dGenesMissing

#############################################################
# parseHOGs: reads the HOG TSV and collects gene counts
#############################################################

def generate_pav_file(hogsfile, outdir, prefix, clade_filter=None):
    """
    Generate a Presence/Absence Variant (PAV) TSV file from HOGs file.
    
    Args:
        hogsfile (str): Path to input HOGs TSV file
        outdir (str): Output directory
        prefix (str): Prefix for output filename
        clade_filter (set, optional): Set of species to include. If None, include all.
    
    Returns:
        str: Path to the generated PAV file
    """
    
    df = pd.read_csv(hogsfile, sep='\t')
    columns_to_drop = ['Gene', 'Tree', 'Parent', 'Clade']
    columns_to_drop = [col for col in columns_to_drop if col in df.columns]
    if columns_to_drop:
        df = df.drop(columns=columns_to_drop)
    
    if clade_filter is not None:
        species_cols = [col for col in df.columns if col not in ['HOG', 'OG']]
        cols_to_keep = ['HOG', 'OG'] + [col for col in species_cols if col in clade_filter]
        df = df[cols_to_keep]
    
    # Convert to presence/absence (1/0)
    species_cols = [col for col in df.columns if col not in ['HOG', 'OG']]
    for col in species_cols:
        df[col] = df[col].apply(lambda x: 1 if pd.notna(x) and str(x).strip() != '' else 0)
    
    pav_file = os.path.join(outdir, f"{prefix}PAV.tsv")
    df.to_csv(pav_file, sep='\t', index=False)
    return pav_file

def generate_count_matrix(hogsfile, outdir, prefix, clade_filter=None):
    """
    Generate a Count Matrix TSV file from HOGs file.
    
    Args:
        hogsfile (str): Path to input HOGs TSV file
        outdir (str): Output directory
        prefix (str): Prefix for output filename
        clade_filter (set, optional): Set of species to include. If None, include all.
    
    Returns:
        str: Path to the generated Count Matrix file
    """
    df = pd.read_csv(hogsfile, sep='\t')
    
    columns_to_drop = ['Gene', 'Tree', 'Parent', 'Clade']
    columns_to_drop = [col for col in columns_to_drop if col in df.columns]
    if columns_to_drop:
        df = df.drop(columns=columns_to_drop)
    
    # Filter species if clade_filter is provided
    if clade_filter is not None:
        species_cols = [col for col in df.columns if col not in ['HOG', 'OG']]
        cols_to_keep = ['HOG', 'OG'] + [col for col in species_cols if col in clade_filter]
        df = df[cols_to_keep]
    
    species_cols = [col for col in df.columns if col not in ['HOG', 'OG']]
    for col in species_cols:
        df[col] = df[col].apply(lambda x: len(str(x).split(', ')) if pd.notna(x) and str(x).strip() != '' else 0)
    
    count_file = os.path.join(outdir, f"{prefix}CountMatrix.tsv")
    df.to_csv(count_file, sep='\t', index=False)
    return count_file

def parseHOGs(hogsfile, clade_filter=None):
    dSpecies = {}
    dGeneNumbers = {}
    dHOGs = {}
    ddHOGs = {}

    with open(hogsfile, "r") as fdh:
        reader = csv.reader(fdh, delimiter="\t")
        header = next(reader)
        for i in range(3, len(header)):
            sp_name = header[i]
            if clade_filter and sp_name not in clade_filter:
                continue
            dSpecies[i] = sp_name

        for line in reader:
            HOGID = line[0]
            OGID = line[1]
            clade_val = line[2]
            gene_counts = []
            filtered_line = []
            for idx in sorted(dSpecies.keys()):
                cell = line[idx] if idx < len(line) else ""
                count = 0 if cell == "" else len(cell.split(", "))
                gene_counts.append(count)
                filtered_line.append(cell)
            dGeneNumbers[HOGID] = gene_counts
            dHOGs[HOGID] = filtered_line

            hog_obj = HOG(HOGID)
            hog_obj.addClade(clade_val)
            for idx in sorted(dSpecies.keys()):
                cell = line[idx] if idx < len(line) else ""
                hog_obj.addMembers(dSpecies[idx], '', cell)
                hog_obj.addCultivar(dSpecies[idx])
            ddHOGs[HOGID] = hog_obj

    return dGeneNumbers, dHOGs, dSpecies, ddHOGs

##############################################
# Extract private (genotype-specific) genes
##############################################

def extract_private_genes(gt_hogs_file, outdir=".", prefix=""):
    with open(gt_hogs_file, "r") as f:
        lines = [x.strip() for x in f if x.strip()]
    if not lines:
        print(f"[INFO] No genotype-specific HOGs found in {gt_hogs_file}.")
        return

    header = lines[0].split("\t")
    species_list = header[1:]
    private_hog_lines = {sp: [] for sp in species_list}
    private_genes = {sp: [] for sp in species_list}

    for line in lines[1:]:
        cols = line.split("\t")
        hog_id = cols[0]
        sp_cols = cols[1:]
        non_empty = [i for i, val in enumerate(sp_cols) if val.strip() != ""]
        if len(non_empty) == 1:
            idx = non_empty[0]
            sp_name = species_list[idx]
            gene_list_str = sp_cols[idx].strip()
            private_hog_lines[sp_name].append(line)
            genes = [g.strip() for g in gene_list_str.split(",")]
            private_genes[sp_name].extend(genes)

    for sp in species_list:
        hog_tsv_name = os.path.join(outdir, f"{prefix}private_HOGs_{sp}.tsv")
        genes_txt_name = os.path.join(outdir, f"{prefix}private_genes_{sp}.txt")
        with open(hog_tsv_name, "w") as fout:
            fout.write("\t".join(header) + "\n")
            for hog_line in private_hog_lines[sp]:
                fout.write(hog_line + "\n")
        with open(genes_txt_name, "w") as fout:
            for g in private_genes[sp]:
                fout.write(g + "\n")
        print(f"[INFO] Species '{sp}': {len(private_hog_lines[sp])} private HOG lines, {len(private_genes[sp])} private genes.")

    summary_file = os.path.join(outdir, f"{prefix}private_gene_counts.txt")
    with open(summary_file, "w") as summary_out:
        summary_out.write("Species\tPrivateGeneCount\n")
        for sp in species_list:
            summary_out.write(f"{sp}\t{len(private_genes[sp])}\n")
    print(f"[INFO] Wrote '{summary_file}' with total private genes per species.")

##################################################
# Build a pan-proteome (merging all categories)
##################################################

def build_pan_proteome(coreHOGs, scHOGs, shellHOGs, gtHOGs,
                       dHOGs, dSpecies, fasta_dir, species_filter, outdir, prefix):
    all_genes_to_write = set()
    cat_dict = {"core": coreHOGs, "single_copy": scHOGs, "shell": shellHOGs, "gt_specific": gtHOGs}
    for hog_dict in cat_dict.values():
        for hog_id, hog_obj in hog_dict.items():
            for sp_name in species_filter:
                if sp_name in hog_obj.members:
                    gene_list_str = hog_obj.members[sp_name].get('', "")
                    if gene_list_str.strip():
                        genes = [g.strip() for g in gene_list_str.split(",")]
                        all_genes_to_write.update(genes)
    proteome_fasta = os.path.join(outdir, f"{prefix}pan_proteome.fa")
    with open(proteome_fasta, "w") as fout:
        for fasta in glob.glob(os.path.join(fasta_dir, "*.fa")) + glob.glob(os.path.join(fasta_dir, "*.fasta")):
            sp = os.path.basename(fasta).split('.')[0]
            if sp not in species_filter:
                continue
            seq_dict = parse_fasta(fasta)
            for gene_id in all_genes_to_write:
                if gene_id in seq_dict:
                    fout.write(f">{gene_id}\n{seq_dict[gene_id]}\n")
    print(f"[INFO] Wrote pan-proteome FASTA with {len(all_genes_to_write)} sequences to: {proteome_fasta}")

##################################################
# Plot gene-count variation heatmap with subplots
##################################################

def plot_genevar_heatmap(dGeneNumbers, dSpecies, coreHOGs, shellHOGs, gtHOGs,
                         outdir, prefix, genevar_filter, transformation):
    sorted_keys = sorted(dSpecies.keys())
    sp_order = [dSpecies[i] for i in sorted_keys]
    hog_ids = sorted(dGeneNumbers.keys())
    data_matrix = [dGeneNumbers[h] for h in hog_ids]
    df = pd.DataFrame(data_matrix, index=hog_ids, columns=sp_order)

    if genevar_filter:
        keep_cols = [c for c in df.columns if c in genevar_filter]
        df = df[keep_cols]

    def is_invariant(row):
        return (row.max() - row.min()) < 2
    df_filtered = df[~df.apply(is_invariant, axis=1)]
    if df_filtered.empty:
        print("[WARNING] No variable data remains for heatmap after filtering.")
        return

    if transformation == "log":
        df_trans = np.log2(df_filtered + 1)
    elif transformation == "zscore":
        df_log = np.log2(df_filtered + 1)
        df_trans = df_log.sub(df_log.mean(axis=1), axis=0).div(df_log.std(axis=1), axis=0).fillna(0)
    else:
        df_trans = df_filtered.copy()

    core_ids = set(coreHOGs.keys())
    shell_ids = set(shellHOGs.keys())
    private_ids = set(gtHOGs.keys())

    core_idx = list(core_ids.intersection(df_trans.index))
    shell_idx = list(shell_ids.intersection(df_trans.index))
    private_idx = list(private_ids.intersection(df_trans.index))

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(18, 6), sharey=False)

    def plot_sub_heatmap(df_subset, title, ax):
        if df_subset.empty:
            ax.text(0.5, 0.5, f"No data for {title}", ha='center', va='center', fontsize=12)
            ax.set_axis_off()
        else:
            sns.heatmap(df_subset, ax=ax, cmap='viridis', cbar_kws={'label': 'Value'})
            ax.set_title(title)
            ax.set_xlabel("Species")
            ax.set_ylabel("HOG ID")

    df_core = df_trans.loc[core_idx] if core_idx else pd.DataFrame()
    df_shell = df_trans.loc[shell_idx] if shell_idx else pd.DataFrame()
    df_private = df_trans.loc[private_idx] if private_idx else pd.DataFrame()

    plot_sub_heatmap(df_core, "Core HOGs", axes[0])
    plot_sub_heatmap(df_shell, "Shell HOGs", axes[1])
    plot_sub_heatmap(df_private, "Private HOGs", axes[2])

    plt.tight_layout()
    for ext in ['png', 'pdf', 'svg']:
        heatmap_file = os.path.join(outdir, f"{prefix}genevar_heatmap.{ext}")
        plt.savefig(heatmap_file, dpi=300)
    plt.close()

    print(f"[INFO] Saved gene variation heatmap to: {prefix}genevar_heatmap.[png|pdf|svg]")

##################################################
# Saturation Analysis (missing from original)
##################################################

def run_saturation_analysis(ddHOGs, species_list, outdir, prefix, bootstrap=100, 
                           marker_core='o', marker_pan='o'):
    """
    Run saturation analysis for the given species list.
    """
    n = len(species_list)
    results = {k: [] for k in range(1, n+1)}
    
    for k in range(1, n+1):
        for _ in range(bootstrap):
            subset = random.sample(species_list, k)
            core_count, pan_count = classify_subset_pan_core(subset, ddHOGs)
            results[k].append((core_count, pan_count))
    
    k_vals, core_means, core_stds, pan_means, pan_stds = [], [], [], [], []
    for k in range(1, n+1):
        core_list = [x[0] for x in results[k]]
        pan_list = [x[1] for x in results[k]]
        k_vals.append(k)
        core_means.append(statistics.mean(core_list))
        core_stds.append(statistics.pstdev(core_list))
        pan_means.append(statistics.mean(pan_list))
        pan_stds.append(statistics.pstdev(pan_list))
    
    plt.figure(figsize=(7, 6))
    plt.errorbar(k_vals, core_means, yerr=core_stds, color='red',
                 label='Core', fmt=f'-{marker_core}')
    plt.errorbar(k_vals, pan_means, yerr=pan_stds, color='blue',
                 label='Shell+Private', fmt=f'-{marker_pan}')
    
    plt.xlabel("Number of Accessions")
    plt.ylabel("Number of HOGs")
    plt.title("Saturation Analysis")
    plt.legend()
    plt.tight_layout()
    
    for ext in ['png', 'pdf', 'svg']:
        sat_file = os.path.join(outdir, f"{prefix}saturation_analysis.{ext}")
        plt.savefig(sat_file, dpi=300)
    plt.close()
    
    print(f"[INFO] Saved saturation analysis figure to: {prefix}saturation_analysis.[png|pdf|svg]")

##################################################
# Saturation Analysis for defined clades (cladepair)
##################################################

def classify_subset_pan_core(subset_species, ddHOGs):
    """
    For a given subset of species, classify each HOG as:
      - core: present in all species in the subset
      - shell+private: present in >=1 but < subset_size
    Returns (core_count, shell_private_count).
    """
    core_count = 0
    shell_private_count = 0
    total = len(subset_species)
    for hog_id, hog_obj in ddHOGs.items():
        presence_count = 0
        for sp in subset_species:
            if sp in hog_obj.members:
                gene_list_str = hog_obj.members[sp].get('', "")
                if gene_list_str.strip():
                    presence_count += 1
        if presence_count == total:
            core_count += 1
        elif 1 <= presence_count < total:
            shell_private_count += 1
    return core_count, shell_private_count

def run_saturation_analysis_defined(ddHOGs, clade1, clade2, outdir, prefix,
                                    bootstrap=100, 
                                    marker_core_clade1='^', marker_pan_clade1='^',
                                    marker_core_clade2='s', marker_pan_clade2='s',
                                    color_core_clade1='#c0392b', color_pan_clade1='#f1c40f',
                                    color_core_clade2='#c0392b', color_pan_clade2='#3498db'):
    # Saturation for clade1
    n1 = len(clade1)
    results1 = {k: [] for k in range(1, n1+1)}
    for k in range(1, n1+1):
        for _ in range(bootstrap):
            subset = random.sample(clade1, k)
            core_count, pan_count = classify_subset_pan_core(subset, ddHOGs)
            results1[k].append((core_count, pan_count))
    k_vals1, core_means1, core_stds1, pan_means1, pan_stds1 = [], [], [], [], []
    for k in range(1, n1+1):
        core_list = [x[0] for x in results1[k]]
        pan_list  = [x[1] for x in results1[k]]
        k_vals1.append(k)
        core_means1.append(statistics.mean(core_list))
        core_stds1.append(statistics.pstdev(core_list))
        pan_means1.append(statistics.mean(pan_list))
        pan_stds1.append(statistics.pstdev(pan_list))
        
    # Saturation for clade2
    n2 = len(clade2)
    results2 = {k: [] for k in range(1, n2+1)}
    for k in range(1, n2+1):
        for _ in range(bootstrap):
            subset = random.sample(clade2, k)
            core_count, pan_count = classify_subset_pan_core(subset, ddHOGs)
            results2[k].append((core_count, pan_count))
    k_vals2, core_means2, core_stds2, pan_means2, pan_stds2 = [], [], [], [], []
    for k in range(1, n2+1):
        core_list = [x[0] for x in results2[k]]
        pan_list  = [x[1] for x in results2[k]]
        k_vals2.append(k)
        core_means2.append(statistics.mean(core_list))
        core_stds2.append(statistics.pstdev(core_list))
        pan_means2.append(statistics.mean(pan_list))
        pan_stds2.append(statistics.pstdev(pan_list))
        
    plt.figure(figsize=(7, 6))
    # For clade1
    plt.errorbar(k_vals1, core_means1, yerr=core_stds1, color=color_core_clade1,
                 label='Clade1 Core', fmt=f'-{marker_core_clade1}')
    plt.errorbar(k_vals1, pan_means1, yerr=pan_stds1, color=color_pan_clade1,
                 label='Clade1 Pan', fmt=f'-{marker_pan_clade1}')
    # For clade2
    plt.errorbar(k_vals2, core_means2, yerr=core_stds2, color=color_core_clade2,
                 label='Clade2 Core', fmt=f'--{marker_core_clade2}')
    plt.errorbar(k_vals2, pan_means2, yerr=pan_stds2, color=color_pan_clade2,
                 label='Clade2 Pan', fmt=f'--{marker_pan_clade2}')

    plt.xlabel("Number of Accessions")
    plt.ylabel("Number of HOGs")
    plt.title("Saturation Analysis by Clade")
    plt.legend()
    plt.tight_layout()
    for ext in ['png', 'pdf', 'svg']:
        sat_file = os.path.join(outdir, f"{prefix}saturation_analysis_defined.{ext}")
        plt.savefig(sat_file, dpi=300)
    plt.close()
    print(f"[INFO] Saved defined saturation analysis figure to: {prefix}saturation_analysis_defined.[png|pdf|svg]")

##################################################
# Main function
##################################################

def main():
    parser = argparse.ArgumentParser(
        description="Phylogeny-Aware Pangenome Classification Toolkit (PanHOG)"
    )
    
    # Configuration file option
    parser.add_argument("--config", type=str, default="config.yaml",
                        help="Path to YAML configuration file (default: config.yaml)")
    
    # Main options
    parser.add_argument("--pan", action="store_true",
                        help="Run overall (global) pangenome classification (default if --clade not provided).")
    parser.add_argument("--clade", type=str, default=None,
                        help="Comma-separated list of species for clade-specific analysis.")
    parser.add_argument("--proteome", type=str, nargs='?', const='ALL', default=None,
                        help="Build a pan-proteome of specified species (comma-separated). If omitted, includes all final species.")
    parser.add_argument("--genevar", type=str, nargs='?', const='ALL', default=None,
                        help="Plot a heatmap of gene variation. If omitted, includes all final species. Otherwise, specify comma-separated species.")
    parser.add_argument("--saturation", action="store_true",
                        help="Perform bootstrapped saturation analysis (core vs shell+private) in incremental steps.")
    parser.add_argument("--saturation-cladepair", action="store_true",
                        help="Perform saturation analysis for two defined clades and plot them in one figure.")
    parser.add_argument("-b", "--bootstrap", type=int, default=100,
                        help="Number of random combinations for saturation analysis (default=100).")
    # Optional markers for global saturation analysis
    parser.add_argument("--marker-core", default='o',
                        help="Marker shape for Core line in saturation analysis (default: 'o').")
    parser.add_argument("--marker-pan", default='o',
                        help="Marker shape for Shell+Private line in saturation analysis (default: 'o').")
    # Options for cladepair saturation analysis
    parser.add_argument("--clade1", type=str,
                        default="Arabis_alpina,ET_AA21_2,ET_AA6,ET_AA7,ET_AA14,ET_AA23,ET_AA22",
                        help="Comma-separated list of species for clade1 (default provided).")
    parser.add_argument("--clade2", type=str,
                        default="Col_PEK,ET108_1,ET131_2,ET133_2,ET148_10,ET33_1,ET49_2,ET53_7,ET3_1,ET96_1,ET105_1,ET173_3,ET150_1",
                        help="Comma-separated list of species for clade2 (default provided).")
    parser.add_argument("--marker-core-clade1", default='^',
                        help="Marker shape for Clade1 Core line (default: '^').")
    parser.add_argument("--marker-pan-clade1", default='^',
                        help="Marker shape for Clade1 Pan line (default: '^').")
    parser.add_argument("--marker-core-clade2", default='s',
                        help="Marker shape for Clade2 Core line (default: 's').")
    parser.add_argument("--marker-pan-clade2", default='s',
                        help="Marker shape for Clade2 Pan line (default: 's').")
    parser.add_argument("--color-core-clade1", default='#c0392b',
                        help="Color for Clade1 Core line (default: '#c0392b').")
    parser.add_argument("--color-pan-clade1", default='#f1c40f',
                        help="Color for Clade1 Pan line (default: '#f1c40f').")
    parser.add_argument("--color-core-clade2", default='#c0392b',
                        help="Color for Clade2 Core line (default: '#c0392b').")
    parser.add_argument("--color-pan-clade2", default='#3498db',
                        help="Color for Clade2 Pan line (default: '#3498db').")
    parser.add_argument("--pav", action="store_true",
                        help="Generate Presence/Absence Variant (PAV) TSV file.")
    parser.add_argument("--matrix", action="store_true",
                        help="Generate Count Matrix TSV file.")
    
    # Input files (required only if not using config file)
    parser.add_argument("--hog", help="Path to HOGs TSV file (e.g. N0.tsv)")
    parser.add_argument("--fasta", help="Directory containing FASTA files")
    # Output options
    parser.add_argument("-o", "--output", type=str, default=".",
                        help="Output directory to write files (default: current directory)")
    parser.add_argument("-p", "--prefix", type=str, default="",
                        help="Prefix to add to output file names (default: none)")
    # Genevar transformation options
    parser.add_argument("--zscore", action="store_true",
                        help="Apply z-score normalization for gene variation heatmap (used with --genevar)")
    parser.add_argument("--log", action="store_true",
                        help="Apply log2(count+1) transformation for gene variation heatmap (used with --genevar)")

    args = parser.parse_args()

    # Load configuration from YAML file
    config = load_config(args.config)
    
    # If config file was loaded, make --hog and --fasta optional
    if config and (args.hog is None or args.fasta is None):
        if 'hog_file' in config and args.hog is None:
            args.hog = config['hog_file']
        if 'fasta_dir' in config and args.fasta is None:
            args.fasta = config['fasta_dir']
    
    if args.hog is None or args.fasta is None:
        parser.error("The following arguments are required: --hog, --fasta (or provide them in config file)")
    
    args = merge_config_with_args(config, args)

    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)

    if not args.pan and not args.clade:
        args.pan = True

    clade_filter = None
    if args.clade:
        clade_filter = set(x.strip() for x in args.clade.split(","))

    hogsfile = args.hog
    fasta_dir = args.fasta
    outdir = args.output
    prefix = args.prefix

    genevar_filter = None
    if args.genevar:
        if args.genevar == 'ALL':
            genevar_filter = None
        else:
            genevar_filter = set(x.strip() for x in args.genevar.split(","))

    proteome_filter = None
    if args.proteome:
        if args.proteome == 'ALL':
            proteome_filter = None
        else:
            proteome_filter = set(x.strip() for x in args.proteome.split(","))

    if args.zscore and args.log:
        sys.exit("[ERROR] Specify only one of --zscore or --log for gene variation heatmap transformation.")
    if args.genevar:
        transformation = "zscore" if args.zscore else "log"
    else:
        transformation = None

    # Generate PAV and/or Count Matrix files if requested
    if args.pav:
        pav_file = generate_pav_file(hogsfile, outdir, prefix, clade_filter)
        print(f"[INFO] Generated PAV file: {pav_file}")
    
    if args.matrix:
        count_file = generate_count_matrix(hogsfile, outdir, prefix, clade_filter)
        print(f"[INFO] Generated Count Matrix file: {count_file}")
    
    # Parse HOGs and missing genes (global or clade-specific)
    dGeneNumbers, dHOGs, dSpecies, ddHOGs = parseHOGs(hogsfile, clade_filter)
    dGenesMissing = getMissingGenes(hogsfile, fasta_dir, clade_filter)
    cloud = []
    for sp in dGenesMissing:
        for g in dGenesMissing[sp]:
            cloud.append(g)

    # Classify HOGs
    coreHOGs = {}
    scHOGs = {}
    shellHOGs = {}
    gtHOGs = {}

    for hog_id, counts in dGeneNumbers.items():
        total_species = len(dSpecies)
        if counts.count(0) == 0:
            coreHOGs[hog_id] = ddHOGs[hog_id]
        if counts.count(1) == total_species:
            scHOGs[hog_id] = ddHOGs[hog_id]
        if counts.count(0) in range(1, total_species - 1) and sum(counts) != 1:
            shellHOGs[hog_id] = ddHOGs[hog_id]
        if counts.count(0) != 0 and sum(counts) == 1:
            gtHOGs[hog_id] = ddHOGs[hog_id]
        if counts.count(0) == total_species - 1:
            gtHOGs[hog_id] = ddHOGs[hog_id]

    print("Core HOGs:", len(coreHOGs))
    print("Single-copy HOGs:", len(scHOGs))
    print("Shell HOGs:", len(shellHOGs))
    print("Genotype-specific HOGs:", len(gtHOGs))
    print("Cloud (unassigned genes):", len(cloud))
    print("Total classified HOGs:", len(coreHOGs) + len(shellHOGs) + len(gtHOGs))

    sorted_species_indices = sorted(dSpecies.keys())
    final_species_list = [dSpecies[i] for i in sorted_species_indices]

    with open(os.path.join(outdir, f"{prefix}core.HOGs.tsv"), 'w') as fout:
        print('hog', *final_species_list, sep='\t', file=fout)
        for hog_id in coreHOGs:
            print(hog_id, *dHOGs[hog_id], sep='\t', file=fout)

    with open(os.path.join(outdir, f"{prefix}single-copy.HOGs.tsv"), 'w') as fout:
        print('hog', *final_species_list, sep='\t', file=fout)
        for hog_id in scHOGs:
            print(hog_id, *dHOGs[hog_id], sep='\t', file=fout)

    with open(os.path.join(outdir, f"{prefix}shell.HOGs.tsv"), 'w') as fout:
        print('hog', *final_species_list, sep='\t', file=fout)
        for hog_id in shellHOGs:
            print(hog_id, *dHOGs[hog_id], sep='\t', file=fout)

    with open(os.path.join(outdir, f"{prefix}gt-specific.HOGs.tsv"), 'w') as fout:
        print('hog', *final_species_list, sep='\t', file=fout)
        for hog_id in gtHOGs:
            print(hog_id, *dHOGs[hog_id], sep='\t', file=fout)

    with open(os.path.join(outdir, f"{prefix}cloud.unassigned_genes.tsv"), 'w') as fout:
        for sp in dGenesMissing:
            print(sp, *dGenesMissing[sp], sep='\t', file=fout)

    extract_private_genes(os.path.join(outdir, f"{prefix}gt-specific.HOGs.tsv"), outdir, prefix)

    if args.proteome is not None:
        if proteome_filter is None:
            species_filter = set(final_species_list)
        else:
            species_filter = proteome_filter.intersection(set(final_species_list))
        if species_filter:
            build_pan_proteome(coreHOGs, scHOGs, shellHOGs, gtHOGs,
                               dHOGs, dSpecies, fasta_dir, species_filter, outdir, prefix)
        else:
            print("[WARNING] --proteome species do not match final dataset. Skipping pan-proteome.")

    if args.genevar is not None:
        if transformation is None:
            transformation = "log"
        if genevar_filter is None:
            gfilt = set(final_species_list)
        else:
            gfilt = genevar_filter.intersection(set(final_species_list))
            if not gfilt:
                print("[WARNING] --genevar species do not match final dataset. Skipping heatmap.")
                return
        plot_genevar_heatmap(dGeneNumbers, dSpecies, coreHOGs, shellHOGs, gtHOGs,
                             outdir, prefix, gfilt, transformation)

    if args.saturation:
        run_saturation_analysis(ddHOGs, final_species_list, outdir, prefix,
                                bootstrap=args.bootstrap,
                                marker_core=args.marker_core,
                                marker_pan=args.marker_pan)

    if args.saturation_cladepair:
        clade1 = [x.strip() for x in args.clade1.split(",")]
        clade2 = [x.strip() for x in args.clade2.split(",")]
        run_saturation_analysis_defined(ddHOGs, clade1, clade2, outdir, prefix,
                                        bootstrap=args.bootstrap,
                                        marker_core_clade1=args.marker_core_clade1,
                                        marker_pan_clade1=args.marker_pan_clade1,
                                        marker_core_clade2=args.marker_core_clade2,
                                        marker_pan_clade2=args.marker_pan_clade2,
                                        color_core_clade1=args.color_core_clade1,
                                        color_pan_clade1=args.color_pan_clade1,
                                        color_core_clade2=args.color_core_clade2,
                                        color_pan_clade2=args.color_pan_clade2)

if __name__ == "__main__":
    main()

