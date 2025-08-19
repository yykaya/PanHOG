# PanHOG

<table>
<tr>
<td width="400">
  <img src="panhog.png" alt="PanHOG logo" width="480"/>
</td>
<td>
A phylogeny-aware toolkit for classifying and annotating Hierarchical Orthologous Groups (HOGs) in pangenomic datasets. It's a flexible command-line toolkit for classifying Orthofinder derived HOGs across multiple genomes in a pangenome-aware, phylogeny-informed context. It supports core/shell/private gene classification, heatmap visualization, pan-proteome generation, and saturation analysis. Now supports both command-line arguments and external configuration files.
</td>
</tr>
</table>
---

## Key Features

* **Global Pangenome Classification (`--pan`)**
* **Clade-Specific Analysis (`--clade species1,species2`)**
* **Pan-Proteome Construction (`--proteome`)**
* **Gene Variation Heatmap (`--genevar`)**
* **Bootstrapped Saturation Analysis (`--saturation`)**
* **Clade-based Saturation Analysis (`--saturation-cladepair`)**
* **Customizable via Config File (`--config config.yaml`)**
* **Pangene Integration for Gene Presence/Absence Analysis**

---

## Documentation
- [Install dependencies](Installation.md)
- [Configuration Guide](README_config.md) - Detailed guide for configuring PanHOG with YAML files
- [Pangene Integration Guide](README_pangene.md) - Instructions for using the pangene integration module


---

## Basic Usage

```bash
python PanHOG.py --hog N0.tsv --fasta ./peptides/ --pan -o results/ -p run1_
```

### With clade-specific analysis:

```bash
python PanHOG.py --hog N0.tsv --fasta ./peptides/ \
  --clade Arabis_alpina,ET_AA21_2,ET_AA6 --pan -o results/ -p cladeA_
```

### With pan-proteome and gene variation heatmap:

```bash
python PanHOG.py --hog N0.tsv --fasta ./peptides/ \
  --proteome ALL --genevar ALL --zscore -o results/ -p viz_
```

---

## ⚙️ Configuration File

You can now use a YAML config file to set advanced parameters like colors, markers, labels, and clade definitions.

### Sample `config.yaml`

```yaml
bootstrap: 10000
marker_core_clade1: "^"
marker_core_clade2: "s"
marker_pan_clade1: "^"
marker_pan_clade2: "s"
color_core_clade1: "#c0392b"
color_pan_clade1: "#f1c40f"
color_core_clade2: "#c0392b"
color_pan_clade2: "#3498db"
clade1:
  - Arabis_alpina
  - AA1
  - AA2
....
clade2:
  - Col_PEK
  - Col-CEN
....
```

### Run using config:

```bash
python PanHOG.py --hog N0.tsv --fasta ./peptides/ --saturation-cladepair --config config.yaml
```

> Any command-line flag will **override** the corresponding config value.

---

## Output Files

* `core.HOGs.tsv`, `shell.HOGs.tsv`, `gt-specific.HOGs.tsv`, `single-copy.HOGs.tsv`
* `cloud.unassigned_genes.tsv`
* `private_genes_<species>.txt`
* `pan_proteome.fa`
* `genevar_heatmap.[png|pdf|svg]`
![Gene Variability Heatmap](genevar_heatmap.png)
* `saturation_analysis.[png|pdf|svg]`
![Saturation Analysis](Saturation_byClade.png)
---

## Recommendations

* Use `--proteome` to extract FASTA of shared pangenes.
* Use `--saturation-cladepair` for insight into core/pan genome expansion across defined clades.
* Use `--genevar` with `--zscore` for population-scale expansions or contractions.

---

## Contact & Citation

This tool is currently in **beta**. For questions, contributions, or citation requests, please contact the developer or include the GitHub link in your reference.

---


Happy pangenomics with **PanHOG**! 🐼
