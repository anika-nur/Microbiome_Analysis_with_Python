# Microbiome Analysis with Python

A computational pipeline for analyzing the gut microbiome of *Jadera haematoloma* (soapberry bugs) and *Boisea trivittata* (boxelder bugs) using 16S rRNA gene amplicon sequencing data. The project covers end-to-end microbiome bioinformatics — from raw data cleaning through alpha/beta diversity, taxonomic composition, differential abundance testing, and environmental variable modeling — entirely in Python.

---

## Table of Contents

- [Background](#background)
- [Repository Structure](#repository-structure)
- [Data Description](#data-description)
- [Analysis Overview](#analysis-overview)
- [Key Findings](#key-findings)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Methods Summary](#methods-summary)
- [Limitations and Caveats](#limitations-and-caveats)
- [License](#license)
- [Contact](#contact)

---

## Background

Insect-associated microbiomes play important roles in host nutrition, immunity, and ecological adaptation. This project investigates the gut microbial communities of two true bug species — primarily *J. haematoloma* (Jhae, n ≈ 161 after QC) with a small number of *B. trivittata* (Btri, n = 2) — collected across multiple cities in Florida in 2019. Samples were sequenced using the 16S rRNA V4 region and taxonomically classified against the SILVA 138 database.

The analysis asks:

1. Do wing morphology (long-wing vs. short-wing), sex, host plant, or geographic location influence gut microbial diversity?
2. What is the core microbiome of *J. haematoloma*?
3. Can environmental variables (temperature, humidity, wind speed) explain variation in microbial community structure?

---

## Repository Structure

```
Microbiome_Analysis_with_Python/
│
├── microbiome_Anika.ipynb              # Main analysis notebook (120 cells)
│
├── bug.microbiome.metadata.tsv         # Sample metadata (species, morph, sex, city, host plant, environmental vars)
├── 16s.v4.AFiles.feature.table.tsv     # ASV feature table (2,680 ASVs × 204 samples)
├── 16s.v4.AFiles.rep.seqs.fna          # Representative sequences for each ASV (FASTA)
├── 16s.v4.AFiles.silva138.taxonomy.tsv # Taxonomic assignments (SILVA 138) with confidence scores
│
├── SampleMitoReads.xlsx                # Supplementary mitochondrial read counts
├── blast_tracking_template.xlsx        # BLAST tracking template
│
└── README.md
```

---

## Data Description

| File | Description | Dimensions |
|------|-------------|------------|
| `bug.microbiome.metadata.tsv` | Per-sample metadata including species, life stage, wing morph (LW/SW), sex, collection city/site, GPS coordinates, host plant, temperature, humidity, wind speed, and DNA quality metrics | ~408 rows × 20 columns |
| `16s.v4.AFiles.feature.table.tsv` | ASV count matrix exported from a BIOM file; rows are ASV hashes, columns are sample IDs | 2,680 ASVs × 204 samples |
| `16s.v4.AFiles.silva138.taxonomy.tsv` | Taxonomic classification of each ASV at domain through species level, with classifier confidence scores | 2,680 rows |
| `16s.v4.AFiles.rep.seqs.fna` | FASTA file of representative sequences for each ASV | 2,681 sequences |
| `SampleMitoReads.xlsx` | Mitochondrial read counts per sample | — |
| `blast_tracking_template.xlsx` | Template for tracking BLAST queries | — |

---

## Analysis Overview

The notebook is organized into 14 sections:

1. **Metadata Cleaning** — Removes QC-type annotation rows, coerces numeric columns, and retains only samples that passed IMR quality control.
2. **Feature Table Cleaning** — Loads the ASV table, identifies and merges duplicate plate entries (e.g., `FJ190108-017_plate1` / `_plate2`), and aligns sample IDs with metadata.
3. **Taxonomy Table Cleaning** — Parses semicolon-delimited SILVA taxonomy strings into separate rank columns (domain through species) and filters out non-bacterial and mitochondrial/chloroplast ASVs.
4. **ID Matching Validation** — Asserts that sample IDs and ASV IDs are consistent across all three tables.
5. **Sample Size Assessment** — Evaluates group balance (Jhae vs. Btri, LW vs. SW, host plants, cities) to determine which comparisons are statistically feasible.
6. **Sequencing Depth Analysis** — Summarizes per-sample read depths and visualizes their distribution.
7. **Alpha Diversity** — Rarefies to a uniform depth, computes Observed ASVs, Shannon, Simpson, and Pielou's evenness, and compares groups using Mann–Whitney U tests (morph, sex, host plant) and Kruskal–Wallis tests (city).
8. **Taxonomic Composition** — Calculates relative abundance at phylum and genus levels and generates stacked bar plots.
9. **Core Microbiome** — Identifies genera present in ≥50% of *J. haematoloma* samples and visualizes prevalence vs. mean abundance.
10. **Host Plant Comparison** — Compares microbial profiles between GRT (golden rain tree) and BV (balloon vine) host plants.
11. **Beta Diversity** — Computes Bray–Curtis distances, performs PCoA ordination, and runs PERMANOVA to test whether morph, sex, host plant, or city structure the community.
12. **Differential Abundance** — Uses Mann–Whitney U tests with FDR correction to identify genera differing between groups.
13. **ANCOM** — Implements Analysis of Composition of Microbiomes from scratch to handle compositional data; applies it to morph and host-plant comparisons.
14. **Environmental Variables** — Correlates temperature, humidity, and wind speed with diversity metrics using Spearman correlations, and fits OLS models to ask whether city effects persist after controlling for climate.

---

## Key Findings

- **Wing morph (LW vs. SW) and sex do not significantly affect alpha or beta diversity** in *J. haematoloma*.
- **Host plant matters**: the small group of bugs collected from balloon vine show higher Shannon diversity and a markedly different community composition compared to golden rain tree bugs (confirmed by PERMANOVA and ANCOM).
- **City of collection significantly structures the microbiome** (PERMANOVA), and this effect persists after accounting for temperature, humidity, and wind speed.
- **The core microbiome is small**: only ~10 genera exceed 50% prevalence, dominated by an unclassified Rhizobiaceae, *Wolbachia*, and *Lactococcus*.
- **Temperature is positively correlated with Shannon diversity**; humidity and wind speed show no significant relationship.

---

## Dependencies

The analysis requires the following Python packages:

- **Python** ≥ 3.8
- [pandas](https://pandas.pydata.org/)
- [NumPy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)
- [SciPy](https://scipy.org/)
- [scikit-bio](http://scikit-bio.org/) — for alpha/beta diversity, PCoA, PERMANOVA, rarefaction
- [statsmodels](https://www.statsmodels.org/) — for OLS regression and multiple-testing correction

---

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/anika-nur/Microbiome_Analysis_with_Python.git
   cd Microbiome_Analysis_with_Python
   ```

2. **Create and activate a virtual environment (recommended):**
   ```bash
   python -m venv myenv
   source myenv/bin/activate        # macOS / Linux
   myenv\Scripts\activate           # Windows
   ```

3. **Install dependencies:**
   ```bash
   pip install pandas numpy matplotlib seaborn scipy scikit-bio statsmodels openpyxl
   ```

4. **Launch Jupyter:**
   ```bash
   pip install jupyter
   jupyter notebook microbiome_Anika.ipynb
   ```

---

## Usage

Open `microbiome_Anika.ipynb` in Jupyter Notebook or JupyterLab and run the cells sequentially. All input data files are expected to be in the same directory as the notebook (the `UPLOADS` variable at the top of the notebook controls this path).

The notebook is self-contained: each section includes markdown explanations of the rationale, the code, and an interpretation of the results. Figures are generated inline.

---

## Methods Summary

| Step | Method | Tool / Function |
|------|--------|-----------------|
| Rarefaction | Random subsampling without replacement to uniform depth | `skbio.stats.subsample_counts` |
| Alpha diversity | Observed ASVs, Shannon, Simpson, Pielou's evenness | `skbio.diversity.alpha_diversity` |
| Alpha group tests | Mann–Whitney U (two groups), Kruskal–Wallis (multiple groups) | `scipy.stats` |
| Beta diversity | Bray–Curtis dissimilarity | `skbio.diversity.beta_diversity` |
| Ordination | Principal Coordinates Analysis (PCoA) | `skbio.stats.ordination.pcoa` |
| Community-level tests | PERMANOVA (999 permutations) | `skbio.stats.distance.permanova` |
| Differential abundance | Mann–Whitney U with Benjamini–Hochberg FDR correction | `scipy.stats`, `statsmodels.stats.multitest` |
| Compositional analysis | ANCOM (custom implementation) | Custom function in notebook |
| Environmental modeling | Spearman correlation, OLS regression | `scipy.stats.spearmanr`, `statsmodels` |

---

## Limitations and Caveats

- **Severe class imbalance**: Jhae (n ≈ 161) vs. Btri (n = 2), and GRT (n ≈ 150+) vs. non-GRT (n = 4). Btri comparisons are descriptive only; host-plant results should be interpreted cautiously.
- **Rarefaction depth**: a minimum depth threshold of 4,500 reads was used; samples below this were excluded.
- **ANCOM implementation**: the notebook includes a from-scratch ANCOM implementation rather than using an established package; results are consistent with Mann–Whitney findings but should be cross-validated for publication.
- **Observational design**: environmental variables (temperature, humidity) are confounded with city and collection date, limiting causal inference.

--- 

## Acknowledgement

This project was conducted under the supervision of Dr. Dave Angelini [Professor of Biology and The Charles C. and Pamela W. Leighton Research Fellow in Biology; Associate Chair of Biology at Colby College, ME, USA]. His guidance shaped the experimental design, analytical approach, and biological interpretation throughout this work.
