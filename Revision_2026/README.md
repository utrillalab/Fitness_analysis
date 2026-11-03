# Revision 2026

Code and data for the revised version of Hidalgo et al. (2026), submitted to *npj Systems Biology and Applications*.
The scripts in this folder correct and extend the analyses of the original notebooks (`../Code_Notebooks/`,
`../Functional_Analysis/`, `../ME_model/`) and produce every figure and Supplementary Data file of the revised article.

Main changes with respect to the original analysis:

- **Fitness values.** The original input files averaged replicates of different conditions, and the glucosamine
  column was taken from N-acetylglucosamine (set1IT061) instead of D-glucosamine hydrochloride (set1IT063), the
  condition of the proteomics of Schmidt et al. (2016). `regenerate_inputs.py` rebuilds the inputs with correctly
  paired replicates; the corrected files replace the originals in `../Functional_Analysis/Infiles/`.
- **GO enrichment.** ViSEAGO was rerun with the corrected inputs (`R1_analyses/a15_viseago.R`); the results replace
  those in `../Functional_Analysis/Viseago_Res/`.
- **New analyses** requested by the reviewers: measurement error (t statistics of Price et al. 2018), proteome mass
  by category, ME-model predictions, LTEE loss-of-function mutations, *P. putida* expression, regulatory
  redundancy, GO enrichment across nine species and structure-based annotation of uncharacterized genes.

## Folder layout

| Path | Content |
|---|---|
| `*.py` (this folder) | Corrected inputs, Table S1 and Figs. 1, 2 and 6 |
| `R1_analyses/` | All other analyses, figures and tables (`aNN_*.py`, `aNN_*.R`); results go to `R1_analyses/output/` |
| `R1_analyses/input_DLT/` | Regulatory network and null-model code of D. Ledezma-Tejeida (Fig. 6, Supplementary Fig. 4) |
| `regenerated_inputs/` | Corrected fitness inputs written by `regenerate_inputs.py` |
| `external/` | Material from co-authors and the original supplementary tables (see below) |
| `original_submission/` | Input files and ViSEAGO results of the first submission (from the tag `original-submission`), kept only for `regenerate_inputs.py` and for the comparisons with the original analysis (`verify_tableS1.py`, `a09`, `a09b`, `a15` setting `original`, `a15b`, `a15c`) |
| `Input_Data_Price2018/`, `Input_Data_Putida/` | Public data, **not included**; download them as described below |

The final figures are in `../Output/` and the Supplementary Data files in `../Supplementary_Data/`.

## Requirements

- Python 3.12 with pandas, numpy, scipy, matplotlib, openpyxl and Pillow (tested with pandas 3.0.6, numpy 2.5.3,
  scipy 1.18.1, matplotlib 3.11.2, openpyxl 3.1.5, Pillow 12.3.0). The figures use the Arial font.
- R 4.6 with ViSEAGO (1.26.0, Bioconductor; it installs topGO 2.64.0), dplyr, readxl and svglite, for the GO
  enrichment (`a15_viseago.R`, `a15c_fig3_dotplot.R`, `a29_table_S14.R`, `a30_table_S15_GO.R`).
- `a36_fig7_npj.py` renders the PNG of Fig. 7 with a headless Chromium browser; set the environment variable `EDGE`
  to its path if Microsoft Edge is not at the default Windows location.
- The ME-model was solved with COBRAme and solvemepy (Python 3.6); see `../ME_model/simulation_conditions.md`.

All paths are relative to the repository. To use the data of another clone, set `FA_REPO` to its root.

## Public input data

These files are public and are not redistributed here. Download them and place them at the paths shown before
running the scripts.

| Path | Source |
|---|---|
| `Input_Data_Price2018/fit_logratios_good.tab`, `fit_t.tab`, `expsUsed` | Price et al. (2018), supplementary data for *E. coli* BW25113 (Keio): https://genomics.lbl.gov/supplemental/bigfit/html/Keio/ (downloaded 6 Oct 2026) |
| `Input_Data_Putida/fit_organism_Putida.tsv`, `t_organism_Putida.tsv`, `exp_organism_Putida.txt`, `organism_Putida_genes.tab` | Fitness Browser, *Pseudomonas putida* KT2440 (orgId `Putida`), downloads page: https://fit.genomics.lbl.gov |
| `external/Pp_FigureS5_share/data/raw/fit_organism_Putida.tsv` | Same file as above |
| `external/Pp_FigureS5_share/data/raw/GCA_000007565_feature_table.txt` | NCBI Assembly GCA_000007565 (*P. putida* KT2440), feature table |
| `external/Pp_FigureS5_share/data/raw/GSE85475_KBojanovic_Processed_data_Stress.xlsx` | GEO GSE85475 (Bojanovic et al. 2017), supplementary file: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85475 |
| `external/essential_lists/41589_2020_593_MOESM2_ESM.xlsx` | Lastiri-Pancardo et al. (2020), Supplementary Data 1 (sheet "Table 1": core and conditional essential genes): https://doi.org/10.1038/s41589-020-0593-y |
| `external/LTEE-Ecoli/` | `git clone https://github.com/barricklab/LTEE-Ecoli.git external/LTEE-Ecoli`, then `git -C external/LTEE-Ecoli checkout 2307710` (curated GenomeDiffs, version of 23 Jun 2026; CC0) |
| `Input_Data_ALEdb/` | ALEdb (https://aledb.org), experiments "LTEE ARA" (ExpID 350) and "LTEE", mutation tables exported as CSV (`Proj_LTEE_Exp_LTEE_ARA_ExpID350_mut.csv`, Ara+1 to Ara+6, 106 clones; `ALEdb_LTEE_AraMinus_mut.csv`, Ara−1 to Ara−6, 96 clones) |

The data of the original analysis (`../Input_Data/`) come from the Fitness Browser (`fit_organism_Keio.tsv`) and from
Schmidt et al. (2016), Supplementary Table S6 (`41587_2016_BFnbt3418_MOESM18_ESM_Mod.xlsx`, modified).

## How to run

From this folder, in this order (each step needs the outputs of the previous ones):

```bash
# 1. corrected inputs and Table S1
python regenerate_inputs.py
python make_tableS1_xlsx.py
python verify_tableS1.py            # optional: diagnosis of the original Table S1
python check_core_lists.py          # optional: checks the 46-condition gene lists

# 2. Figs. 1, 2 and 6
python fig1cd.py                    # uses fig1c.py
python fig2.py
python fig6.py

# 3. analyses (R1_analyses/)
cd R1_analyses
python a01_distribution_thresholds.py
python a02_mass_by_category.py
python a03_me_specificity_threshold.py
python a05_t_statistic.py           # before a04
python a04_expected_cost_resolution.py
python a07_tableS19_core_column.py  # Supplementary Data 19, needed by a08
python a08_claims_vs_core.py
python a12_tablesS17_S18.py
python a32_tablesS17_S18_v3.py      # Supplementary Data 17-18 (after a12)
python a13_fig_missing_fitness.py
python a14_essential_nofit_list.py
Rscript a15_viseago.R               # GO enrichment per condition (slow)
python a15b_compare_viseago.py
python a15d_tables_S2_S12.py
python a25_table1_S13.py
Rscript a29_table_S14.R
python a29_table_S14.py
python a29_fig5.py
Rscript a30_table_S15_GO.R
python a06_tableS15_core_column.py
python a09_ltee.py plus; python a09_ltee.py minus; python a09_ltee.py all
python a09b_deletion_audit.py
python a09c_match_barrick.py
python a10_fig_S4.py
python a17_go_mass_null.py
python a18_fig_S1_mass_quintiles.py
python a19_fig_regression.py
python a20_fig4.py
python a21_putida_reproduce.py
python a22_fig_S5.py
python a23_fig_S1_revised.py
python a24_fig_S_split.py
python a26_fig_S_me_model.py
python a27_threshold_core.py
python a28_fig_S6_summary.py
python a31_tf_redundancy_null.py
python a33_fig_S4_redundancy.py
python a16_table_S20.py             # collects the statistical tests (Supplementary Data 20)
Rscript a15c_fig3_dotplot.R
python a35_fig3_npj.py
python a36_fig7_npj.py
cd ..

# 4. Supplementary Data 1-25
python build_supp_tables.py
```

## Figures

| Figure | Script | Output |
|---|---|---|
| Fig. 1 | `fig1cd.py` (with `fig1c.py`) | `R1_analyses/output/Fig01_revised.*` |
| Fig. 2 | `fig2.py` | `R1_analyses/output/Fig02_glucosamineHCl.*` |
| Fig. 3 | `a15_viseago.R` → `a15c_fig3_dotplot.R` → `a35_fig3_npj.py` | `R1_analyses/output/Fig03_dotplot_corrected.*` |
| Fig. 4 | `a20_fig4.py` | `R1_analyses/output/Fig04_revised.*` |
| Fig. 5 | `a29_table_S14.R` → `a29_table_S14.py` → `a29_fig5.py` | `R1_analyses/output/Fig05_revised.*` |
| Fig. 6 | `fig6.py` | `R1_analyses/output/Fig06_revised.*` |
| Fig. 7 | `a36_fig7_npj.py`, from the figure of J. A. Farías (`external/Figure7_JAF.svg`) | `R1_analyses/output/Fig07_npj.*` |
| Supplementary Fig. 1 | `a24_fig_S_split.py` | `Fig_S_classification.*` |
| Supplementary Fig. 2 | `a24_fig_S_split.py` (with `a18`, `a23`) | `Fig_S_mass.*` |
| Supplementary Fig. 3 | `a26_fig_S_me_model.py` | `Fig_S_me_model.*` |
| Supplementary Fig. 4 | `a33_fig_S4_redundancy.py` | `Fig_S4_redundancy.*` |
| Supplementary Fig. 5 | `a09_ltee.py` (three runs) → `a10_fig_S4.py` | `Fig_S4_v1.*` |
| Supplementary Fig. 6 | `a22_fig_S5.py` (reproduction of O. A. Aguilar-Vera's notebook in `a21_putida_reproduce.py`) | `Fig_S5_revised.*` |
| Supplementary Fig. 7 | a: notebook of L. Soto-Avila (`external/Liz_GO_nine_species/`); b: `a28_fig_S6_summary.py` | `Fig_S6_summary.*` |
| Supplementary Fig. 8 | Figure of L. Soto-Avila (not generated here) | — |

The output names keep the provisional numbers used during the revision; the table gives the final numbers.

## Tables and Supplementary Data

| Table | Script |
|---|---|
| Table 1 | `a25_table1_S13.py` |
| Supplementary Data 1 (fitness and proteomics per gene) | `regenerate_inputs.py` → `make_tableS1_xlsx.py` |
| Supplementary Data 2–12 (GO enrichment per condition) | `a15_viseago.R` → `a15d_tables_S2_S12.py` |
| Supplementary Data 13 (GO terms in glucose and their proteome mass) | `a25_table1_S13.py` |
| Supplementary Data 14 (GO enrichment of the core mean-effect genes) | `a29_table_S14.R` → `a29_table_S14.py` |
| Supplementary Data 15 (GO enrichment of the 108 fitness-enhancing candidates) | `a30_table_S15_GO.R` → `a06_tableS15_core_column.py` |
| Supplementary Data 16 | Unchanged from the original article (`external/Supp_Tables_1-19_Hidalgo_et_al.xlsx`) |
| Supplementary Data 17–18 (structure-based annotation) | `a12_tablesS17_S18.py` → `a32_tablesS17_S18_v3.py`, from the table of J. A. Farías |
| Supplementary Data 19 (conditions per category for each gene) | `a07_tableS19_core_column.py` |
| Supplementary Data 20 (statistical tests of the revision) | `a16_table_S20.py`, from the outputs of `a01`–`a33` |
| Supplementary Data 21 (genes without fitness data annotated as essential) | `a14_essential_nofit_list.py` → `a15d_tables_S2_S12.py` |
| Supplementary Data 22 (GO terms of the nine-species analysis by theme) | `a28_fig_S6_summary.py` |
| Supplementary Data 23–24 (*P. putida* highly expressed mean-effect genes; COG enrichment) | O. A. Aguilar-Vera's package (`external/Pp_FigureS5_share/results/`) |
| Supplementary Data 25 (genes with loss-of-function mutations in the LTEE) | `a09_ltee.py all` |

`build_supp_tables.py` collects all of them into one workbook and writes one file per table to
`Supplementary_Data/`. Scripts that only test statements of the response to the reviewers (`a01`–`a05`, `a08`,
`a13`, `a17`, `a19`, `a27`, `a31`) write their results to `R1_analyses/output/`, and most of them are collected in
Supplementary Data 20.

## LTEE analysis (Supplementary Fig. 5, Supplementary Data 25)

Loss-of-function (LoF) mutations in the 12 LTEE populations compared with the glucose fitness categories
(`a09_ltee.py`, `a09b_deletion_audit.py`, `a09c_match_barrick.py`, `a10_fig_S4.py`; each `a09_ltee.py` run takes a few
minutes because of the 10,000 permutations per subset).

- **LoF:** nonsense SNPs; coding indels whose length is not a multiple of 3; IS insertions inside a coding region;
  deletions that remove whole or partial genes. Missense, synonymous, in-frame and intergenic mutations are not LoF.
- **Amplification artefacts:** in the curated GenomeDiffs, eight deletion rows carry `within=<AMP>:<copy>`, i.e. they
  lie inside one copy of an amplification and no gene is lost. ALEdb drops this attribute, so these rows are removed
  (`a09b`, `a09c`). They include the only apparent LoF of *rpoS*.
- **Universe:** genes with glucose fitness that are in the EcoCyc K-12 gene table. REL606 gene names are matched to
  K-12 b-numbers by name; IS elements and `ECB_` loci are not mapped.
- **Mutators:** ALEdb only includes clones without elevated point-mutation rates. IS-mutator clones (Ara+1 from
  5,000 generations; one Ara−5 clone at 30,000) are taken from the GenomeDiff metadata and analysed separately.
- **Tests:** hypergeometric test per fitness category; event-level permutation (10,000 random gene sets, single-gene
  events drawn in proportion to gene length, multi-gene deletions as blocks of the same number of consecutive genes);
  Mann-Whitney test of glucose fitness. Subsets: all LoF, point-like LoF, parallel (≥ 2 populations), persistent
  (≥ 2 time points) and without IS-mutator clones.

## Material from co-authors (`external/` and `R1_analyses/input_DLT/`)

| Path | Content | Author |
|---|---|---|
| `external/Supp_Tables_1-19_Hidalgo_et_al.xlsx` | Supplementary Tables of the original submission | Authors |
| `external/essential_lists/` | Gene lists of PEC, DEG and the 276 Keio essential genes with RegulonDB 12 names, and the original script `Fig_suppl01_intersection_essential_genes.py` (it uses server paths; `a14_essential_nofit_list.py` is the runnable version) | L. Soto-Avila |
| `external/Liz_GO_nine_species/` | GO enrichment across nine species (notebooks `EnrichmentGO_several_orgs_Suppl_S6.ipynb`, `categorize_fitness data.ipynb`; results `all_gos.csv`, `filter_imp.csv`) | L. Soto-Avila |
| `external/Pp_FigureS5_share/` | *P. putida* fitness vs expression, with its own README and conda environment | O. A. Aguilar-Vera |
| `R1_analyses/input_DLT/` | RegulonDB 11.2 TF–gene network (`NetWorkTFGene.txt`), network without global TFs (`TRN_noglobals.txt`) and null-model code | D. Ledezma-Tejeida |
| `external/JoseAFarias_fig7/`, `external/Figure7_JAF.svg` | Structure-based functional annotation (v3) and Fig. 7 | J. A. Farías |

## References

- Bojanovic K, et al. (2017). Global transcriptional responses to osmotic, oxidative, and imipenem stress
  conditions in *Pseudomonas putida*. *Appl. Environ. Microbiol.* 83, e03236-16.
- Lastiri-Pancardo G, et al. (2020). A quantitative method for proteome reallocation using minimal regulatory
  interventions. *Nat. Chem. Biol.* 16, 1026–1033.
- Phaneuf PV, et al. (2019). ALEdb 1.0: a database of mutations from adaptive laboratory evolution experimentation.
  *Nucleic Acids Res.* 47, D1164–D1171.
- Price MN, et al. (2018). Mutant phenotypes for thousands of bacterial genes of unknown function. *Nature* 557,
  503–509.
- Schmidt A, et al. (2016). The quantitative and condition-dependent *Escherichia coli* proteome. *Nat. Biotechnol.*
  34, 104–110.
- Tenaillon O, et al. (2016). Tempo and mode of genome evolution in a 50,000-generation experiment. *Nature* 536,
  165–170.
