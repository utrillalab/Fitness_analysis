# Fitness_analysis
Code and data to analyze and generate figures for Hidalgo et al. 2026

## Folders

- `Code_Notebooks/`, `Functional_Analysis/`, `ME_model/`: notebooks of the original analysis. The input files in
  `Functional_Analysis/Infiles/` and the GO results in `Functional_Analysis/Viseago_Res/` were replaced by the corrected
  versions of the revision, and so were `Functional_Analysis/BP_Results_genes_noeffect.xlsx` and
  `BP_Results_genes_enhance.xlsx` (GO enrichment of the core mean-effect genes and of the 108 fitness-enhancing
  candidates; tab-separated text, as written by ViSEAGO; `Revision_2026/R1_analyses/a29_table_S14.R` and
  `a30_table_S15_GO.R`). The ME-model simulation conditions are described in `ME_model/simulation_conditions.md`.
- `Revision_2026/`: scripts and data of the revised version (corrected inputs, new analyses, LTEE loss-of-function
  mutations, all figures and Supplementary Data); see its README.
- `Output/`: final figures. `Supplementary_Data/`: Supplementary Data 1-25.
