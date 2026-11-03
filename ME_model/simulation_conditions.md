# ME-model simulation conditions

This folder contains the genome-scale model of metabolism and gene expression (ME-model) used in the article, the
notebooks that solve it for each carbon source and extract the translation fluxes, and the saved solutions.

## Model

- **File:** `iJL1678b_ML_KO_keff_FbFP.pickle`: the *E. coli* K-12 MG1655 ME-model iJL1678b-ME
  (Liu et al., 2014; built with ECOLIme, https://github.com/SBRG/ECOLIme) in which the effective rate constants
  (keffs) were replaced by the values predicted by machine learning in Heckmann et al. (2018), using the quantitative
  proteomics of Schmidt et al. (2016).
- **Framework:** COBRAme (Lloyd et al., 2018; https://github.com/SBRG/cobrame), Python 3.6.13.
- **Solver:** qMINOS through solvemepy (Ma et al., 2017; https://github.com/SBRG/solvemepy).

## Conditions

| Parameter | Value |
|---|---|
| Objective | maximal growth rate, found by bisection on μ between 0 and 1.2 h⁻¹ (precision 10⁻⁶) |
| Unmodeled protein fraction | 0.36 |
| Growth-associated maintenance (GAM) | 34.98 mmol ATP gDW⁻¹ |
| Non-growth-associated maintenance (nGAM) | 1 mmol ATP gDW⁻¹ h⁻¹ |
| Medium | aerobic minimal medium: the glucose uptake reaction is closed and the uptake reaction of the selected carbon source is opened with bounds of −1000 and 1000 mmol gDW⁻¹ h⁻¹, so that growth is limited by the proteome allocation constraints of the model rather than by substrate uptake. All other exchange reactions keep the default bounds of the model (uptake of oxygen, ammonium, phosphate, sulfate and trace ions allowed). |
| Output | flux of the translation reaction of each of the 1,568 genes in the model, used as its predicted use |

LB medium was not simulated because its composition is undefined.

## Carbon sources and predicted growth rates

| Carbon source | Uptake reaction opened | Predicted growth rate (h⁻¹) | Saved solution (`ME_Output/`) |
|---|---|---|---|
| Glucose | `EX_glc__D_e` | 1.040 | `EX_glc__D_e_34_1_036.pkl` |
| Glycerol | `EX_glyc_e` | 1.038 | `EX_glyc_e_34_1_036.pkl` |
| Fructose | `EX_fru_e` | 1.036 | `EX_fru_e_34_1_036.pkl` |
| Mannose | `EX_man_e` | 1.003 | `EX_man_e_34_1_036.pkl` |
| Galactose | `EX_gal_e` | 0.995 | `EX_gal_e_34_1_036.pkl` |
| Succinate | `EX_succ_e` | 0.993 | `EX_succ_e_34_1_036.pkl` |
| Pyruvate | `EX_pyr_e` | 0.989 | `EX_pyr_e_34_1_036.pkl` |
| Glucosamine | `EX_gam_e` | 0.981 | `EX_gam_e_34_1_036.pkl` |
| Xylose | `EX_xyl__D_e` | 0.852 | `EX_xyl__D_e_34_1_036.pkl` |
| Acetate | `EX_ac_e` | 0.713 | `EX_ac_e_34_1_036.pkl` |

The growth rates are the `biomass_dilution` fluxes of the saved solutions. The file names encode GAM (34), nGAM (1)
and the unmodeled protein fraction (0.36).

## How to run

1. Install COBRAme and solvemepy (the Docker image `sbrg/cobrame` from https://hub.docker.com/r/sbrg/cobrame/ includes
   both).
2. Run `Solving_ME_carbon_sources.ipynb` from this folder. It loads `iJL1678b_ML_KO_keff_FbFP.pickle`, solves the model
   for each carbon source and writes the solutions to `ME_Output/`.
3. Run `ME_Fluxes.ipynb` to extract the translation flux of each gene (`ME_Output/Updated_ME_fluxes.csv`), which is
   used for Fig. 1d and Supplementary Fig. 3.

## References

- Heckmann D, et al. (2018). Machine learning applied to enzyme turnover numbers reveals protein structural correlates
  and improves metabolic models. *Nat. Commun.* 9, 5252.
- Liu JK, et al. (2014). Reconstruction and modeling protein translocation and compartmentalization in *Escherichia
  coli* at the genome-scale. *BMC Syst. Biol.* 8, 110.
- Lloyd CJ, et al. (2018). COBRAme: A computational framework for genome-scale models of metabolism and gene
  expression. *PLoS Comput. Biol.* 14, e1006302.
- Ma D, et al. (2017). Reliable and efficient solution of genome-scale models of Metabolism and macromolecular
  Expression. *Sci. Rep.* 7, 40863.
- Schmidt A, et al. (2016). The quantitative and condition-dependent *Escherichia coli* proteome. *Nat. Biotechnol.*
  34, 104–110.
