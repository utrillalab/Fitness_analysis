# Supplementary Figure S5: fitness vs expression in *Pseudomonas putida* KT2440

Code and data to reproduce Supplementary Figure S5: gene fitness on glucose versus expression,
and the COG categories of highly expressed genes in the mean-effect fitness category.

## Contents

```
.
├── README.md
├── environment.yml                         # conda environment (exact versions)
├── P.putida_Fitness_ExpressionTPM.ipynb    # main analysis (Steps 1-8)
├── scripts/
│   └── tpm_from_counts.py                  # raw counts -> TPM table
├── data/
│   ├── raw/                                # input data, as downloaded
│   └── processed/
│       └── GSE85475_TPM.tsv                # output of scripts/tpm_from_counts.py
└── results/
    ├── FigureS5.pdf                        # 2 pages: panels A-B and gene table
    ├── FigureS5_tableC.tsv                 # gene table (full names)
    ├── FigureS5_A.png, FigureS5_B.png      # panel images used in the PDF
    └── TableSX_COG_enrichment.tsv          # COG enrichment, top 5% vs remaining mean-effect genes
```

## Input data (`data/raw/`)

| File | Content | Source |
|---|---|---|
| `GSE85475_KBojanovic_Processed_data_Stress.xlsx` | Processed RNA-seq data (raw counts, RPKM, CDS coordinates); 21 samples | GEO [GSE85475](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85475), supplementary file |
| `fit_organism_Putida.tsv` | Gene fitness for all *P. putida* experiments | Fitness Browser (https://fit.genomics.lbl.gov) |
| `GCA_000007565_feature_table.txt` | Genome annotation of *P. putida* KT2440 | NCBI Assembly GCA_000007565 |
| `COGs.csv` | COG assignment per locus tag | [to be completed] |

## Method summary

1. **TPM** (`scripts/tpm_from_counts.py`): raw counts of all 21 GSE85475 samples are divided by CDS length
   (|Translation Stop − Translation Start| + 1, from the GEO file) and normalized to 1e6 per sample.
   Genes: `PP_` locus tags, excluding tRNA and rRNA (NCBI annotation).
2. **Expression**: mean TPM of the three control replicates (`C_1_0`, `C_2_0`, `C_3_0`), log10(TPM + 1).
3. **Fitness**: experiment `set6IT057 D-Glucose (C)`. Categories: essential (≤ −2), important (−2, −0.5),
   mean-effect [−0.5, 0.5], fitness-enhancing (> 0.5).
4. **Highly expressed**: top 5% of log10(TPM + 1) among genes with fitness and expression (4,477 genes),
   i.e. log10(TPM + 1) ≥ 2.40 (224 genes; 168 in the mean-effect category).
5. **COG**: genes with a single COG letter are plotted in panel B; the others are listed in the table.
6. **Expression filter** (Step 8): COG composition of all mean-effect genes (no filter), and Fisher's exact test
   per COG group and letter, top 5% vs remaining mean-effect genes, Benjamini–Hochberg FDR within each level.

## How to reproduce

```bash
conda env create -f environment.yml
conda activate pp_figureS5
plotly_get_chrome -y                          # Chrome for kaleido (static image export)
python -m ipykernel install --user --name pp_figureS5

# from this directory
python scripts/tpm_from_counts.py             # -> data/processed/GSE85475_TPM.tsv
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=pp_figureS5 P.putida_Fitness_ExpressionTPM.ipynb
```

All paths are relative: run the commands from the repository root.
