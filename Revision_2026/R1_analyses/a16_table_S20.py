"""Supplementary Table S20: all statistical tests added in the revision (promised in the response to R1-G1).

Collects the result tables of a01-a15 (output/) into one workbook: an index sheet plus one sheet per analysis, each with
the reviewer comment it answers, the test and the source file. Values are copied as computed; nothing is re-analysed,
except the regulator enrichment of the core set chosen on 7 Oct (option B), computed here as in a03_m10.
Output: ../Table_S20_statistical_tests.xlsx
"""
import os
import numpy as np, pandas as pd
from scipy import stats
from common import *


def rd(name, **kw):
    return pd.read_csv(os.path.join(OUT, name), **kw)


# regulator enrichment of the option-B core (fitness > 0.5 in >= 10 of 46 conditions and t > 4 in >= 1), as in a03_m10
fa, names46 = load_fitness_all()
net = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'NetworkRegulatorGene.tsv'), sep='\t', comment='#', header=None)
lower = {n.lower(): n for n in names46.dropna().unique()}
regs = {lower[r.lower()] for r in set(net[1].dropna().astype(str)) if r.lower() in lower}
reg_b = set(names46.index[names46.isin(regs)])
N, K = len(fa), len(reg_b & set(fa.index))
core = rd('a05_enhancing_108_candidates_32_core.csv')
core_b = set(core.loc[core.significant_core, 'bnumber'])
x = len(core_b & reg_b); n = len(core_b)
reg = rd('a05_enhancing_regulators_t.csv')
reg = pd.concat([reg, pd.DataFrame([{'criterion': 'fitness > 0.5 in >= 10 and t > 4 in >= 1 (core set used in the paper)',
                                     'min_conditions': 10, 'n_genes': n, 'n_regulators': x,
                                     'fold_enrichment': round((x / n) / (K / N), 4),
                                     'p_hypergeom': stats.hypergeom.sf(x - 1, N, K, n)}])], ignore_index=True)
# R2-m2: without rpsG and rpsR, the only candidates essential in Keio, PEC and DEG (probable Tn-seq artefacts)
art_b = set(core.loc[core.gene.isin(['rpsG', 'rpsR']), 'bnumber'])
for lab, s in [('fitness > 0.5 in >= 10, without rpsG and rpsR', set(core.bnumber) - art_b),
               ('fitness > 0.5 in >= 10 and t > 4 in >= 1 (core), without rpsG and rpsR', core_b - art_b)]:
    xs = len(s & reg_b)
    reg.loc[len(reg)] = {'criterion': lab, 'min_conditions': 10, 'n_genes': len(s), 'n_regulators': xs,
                         'fold_enrichment': round((xs / len(s)) / (K / N), 4), 'p_hypergeom': stats.hypergeom.sf(xs - 1, N, K, len(s))}

det = rd('a04_detectability_exact_generations.csv')
det = det[det.phi_max == 0.48]

ltee = []
for lab, f in [('12 populations', 'a09_all_ltee_lof_tests.csv'), ('Ara+ populations', 'a09_ltee_lof_tests.csv'),
               ('Ara- populations', 'a09_minus_ltee_lof_tests.csv')]:
    t = rd(f); t.insert(0, 'populations', lab); ltee.append(t)
ltee = pd.concat(ltee, ignore_index=True)
keep = ['populations', 'set', 'categories', 'events', 'genes_in_universe'] + \
       [c for c in ltee.columns if c.startswith(('mean-effect_', 'essential_')) or c.startswith('p_')]
ltee = ltee[[c for c in keep if c in ltee.columns]]

go = rd('a15b_essential_terms_with_vs_without_287.csv')
go = (go.assign(with_287=go.p_with_287.notna(), without_287=go.p_without_287.notna())
        .groupby(['GO.ID', 'term']).agg(conditions_with_287=('with_287', 'sum'), conditions_without_287=('without_287', 'sum'))
        .reset_index().sort_values(['conditions_with_287', 'conditions_without_287'], ascending=False))

S = [  # (sheet, reviewer comment, test, description, table, source)
 ('S20a_fitness_distribution', 'R1-M2', 'Descriptive statistics',
  'Fitness distribution per condition: mode, median, mean, skewness, robust SD of the central peak and % of genes per category.',
  rd('a01_distribution_summary.csv'), 'a01_distribution_summary.csv'),
 ('S20b_threshold_sensitivity', 'R1-M2, R2-M2a', 'Sensitivity analysis',
  'Proteome mass fraction per category (median across conditions) for mean-effect bounds of ±0.3, ±0.5, ±0.7 and ±1.0 and essential cut-offs of -1.5, -2 and -3.',
  rd('a01_threshold_sensitivity_summary.csv', header=[0, 1], index_col=[0, 1]).reset_index(), 'a01_threshold_sensitivity_summary.csv'),
 ('S20c_genes_without_fitness', 'R1-M1, R2-M1', 'Scenario analysis',
  'Proteome mass fraction per category when the measured proteins without fitness data are added to the denominator (Supplementary Fig. S1X).',
  rd('a13_missing_fitness_fractions.csv'), 'a13_missing_fitness_fractions.csv'),
 ('S20d_mass_permutation_null', 'R1-m4', 'Permutation test (10,000 random gene sets of equal size), BH correction',
  'Observed proteome mass fraction of each category vs. random gene sets of the same size.',
  rd('a02_m4_permutation_null.csv'), 'a02_m4_permutation_null.csv'),
 ('S20e_mass_effect_sizes', 'R1-M4b, G1', "Cliff's delta, Cohen's d, Kruskal-Wallis, Holm correction",
  'Effect size of log10 protein mass for each category vs. mean-effect genes.',
  rd('a02_M4b_effect_sizes.csv'), 'a02_M4b_effect_sizes.csv'),
 ('S20f_mass_quintiles', 'R1-M4c', 'Hypergeometric test per mass quintile, BH correction',
  'Over- and under-representation of each category in each protein-mass quintile, per condition (pooled summary in sheet S20f2).',
  rd('a02_M4c_mass_quintile_enrichment.csv'), 'a02_M4c_mass_quintile_enrichment.csv'),
 ('S20f2_mass_quintiles_pooled', 'R1-M4c', 'Summary of S20f', 'Number of conditions with significant over/under-representation.',
  rd('a02_M4c_pooled.csv'), 'a02_M4c_pooled.csv'),
 ('S20g_fitness_mass_regression', 'R1-M4a, R1-M4d, G2', 'Linear regression with condition as covariate, gene-clustered SE',
  'Slope of fitness per 10-fold increase in protein mass, by category (top) and interaction model (bottom rows).',
  pd.concat([rd('a02_M4d_regression_by_category.csv'), rd('a02_M4d_interaction_model.csv')], ignore_index=True),
  'a02_M4d_regression_by_category.csv; a02_M4d_interaction_model.csv'),
 ('S20h_enhancing_vs_mean_effect', 'R1-m9', 'Mann-Whitney (location); Brown-Forsythe and Ansari-Bradley (spread)',
  'Protein mass of fitness-enhancing vs. mean-effect genes per condition.',
  rd('a02_m9_enhancing_vs_mean_effect.csv'), 'a02_m9_enhancing_vs_mean_effect.csv'),
 ('S20i_ME_model_status', 'R1-M3', 'Contingency (odds ratio, Fisher)',
  'ME-model status (absent, zero flux, non-zero flux) of the genes of each fitness category, per condition.',
  rd('a03_M3_me_status_by_category.csv'), 'a03_M3_me_status_by_category.csv'),
 ('S20j_ME_flux_correlations', 'R1-m3', 'Spearman correlation',
  'Predicted translation flux vs. measured protein mass (non-zero-flux genes) and vs. fitness.',
  rd('a03_m3_flux_correlations.csv'), 'a03_m3_flux_correlations.csv'),
 ('S20k_detectability', 'R1-M3', 'Linear proteome-cost model (phi_c = 0.48, Scott et al. 2010)',
  'Protein mass fraction needed for a single gene deletion to reach fitness +0.5, given the number of generations of each experiment.',
  det, 'a04_detectability_exact_generations.csv'),
 ('S20l_condition_specificity', 'R1-M8', 'Jaccard index vs. size-matched null; Mann-Whitney',
  'Overlap of each category between conditions and number of conditions per gene.',
  rd('a03_M8_condition_specificity.csv'), 'a03_M8_condition_specificity.csv'),
 ('S20m_t_reclassification', 'R2-M2a, R1-M2', '|fitness| > 0.5 and |t| > 4 (Price et al. 2018)',
  'Genes per category before and after requiring |t| > 4.',
  rd('a05_reclassification_counts.csv'), 'a05_reclassification_counts.csv'),
 ('S20n_t_proteome_fractions', 'R2-M2a', '|fitness| > 0.5 and |t| > 4',
  'Proteome mass fraction per category with cut-offs only and with the t criterion.',
  rd('a05_proteome_fractions_t.csv', header=[0, 1], index_col=0).reset_index(), 'a05_proteome_fractions_t.csv'),
 ('S20o_t_effect_sizes', 'R2-M2a', "Cliff's delta, Holm correction", 'Effect sizes of log10 mass vs. mean-effect with t-aware categories.',
  rd('a05_effect_sizes_t.csv'), 'a05_effect_sizes_t.csv'),
 ('S20p_enhancing_threshold', 'R1-m10, R2-m1b, R2-m2', 'Hypergeometric test (RegulonDB regulators)',
  'Regulator enrichment for recurrence thresholds of 3-25 of 46 conditions, for the candidates (fitness > 0.5 in >= k '
  'conditions) and for the core (also t > 4 in >= 1 of them), with and without rpsG and rpsR (probable Tn-seq artefacts).',
  rd('a27_threshold_core_sensitivity.csv'), 'a27_threshold_core_sensitivity.csv'),
 ('S20q_enhancing_core', 'R2-M2a, R2-m1b, R2-m2', 'Hypergeometric test (RegulonDB regulators)',
  'Regulator enrichment with the t criterion; row 7 = core set used in the paper (71 genes); last two rows = candidates and core without rpsG and rpsR (essential in Keio, PEC and DEG; probable Tn-seq artefacts).',
  reg, 'a05_enhancing_regulators_t.csv + this script'),
 ('S20r_LTEE_LoF', 'R1-M6, R2-M3', 'Hypergeometric test; event-level permutation (10,000; gene length and deletion blocks preserved); Mann-Whitney',
  'Fitness categories of genes with loss-of-function mutations in the LTEE vs. all genes, for all, point-like, parallel, persistent and non-IS-mutator subsets.',
  ltee, 'a09_all/a09/a09_minus_ltee_lof_tests.csv'),
 ('S20s_GO_without_287', 'R1-M1', 'topGO elim + Fisher (p < 0.01)',
  'Essential-category GO terms with and without the 287 essential genes lacking fitness data: number of conditions (of 11) in which each term is enriched.',
  go, 'a15b_essential_terms_with_vs_without_287.csv'),
 ('S20t_GO_reproducibility', 'R1-M1, R2-M1', 'Jaccard index of enriched GO terms',
  'Overlap of enriched GO terms between the original and regenerated analyses (repo vs. rerun; original vs. corrected fitness; with vs. without the 287 genes).',
  rd('a15b_viseago_comparison.csv'), 'a15b_viseago_comparison.csv'),
 ('S20u_GO_term_mass_null', 'R1-M5b', 'Permutation test (10,000 random sets of the same number of measured proteins; all proteins and same category), BH correction',
  'Observed vs. expected proteome mass of every enriched GO term (corrected ViSEAGO results), per condition and fitness category.',
  rd('a17_go_mass_null.csv'), 'a17_go_mass_null.csv'),
]

# R2-M2b: D. Ledezma-Tejeida's null models for TF redundancy (her Supplementary Tables 1 and 2, input_DLT/tablas_supp.xlsx)
_dlt = pd.read_excel(os.path.join(HERE, 'input_DLT', 'tablas_supp.xlsx'), header=None)
_starts = [i for i, v in enumerate(_dlt[0]) if v in ('Table S1', 'Table S2')]


def _dlt_table(i0):
    t = _dlt.iloc[i0 + 1:].reset_index(drop=True)
    end = next((k for k, v in enumerate(t[0]) if pd.isna(v)), len(t))
    t = t.iloc[:end].dropna(axis=1, how='all')
    t.columns = t.iloc[0]
    return t.iloc[1:].reset_index(drop=True)


S += [
 ('S20v_TF_redundancy_random_sets', 'R2-M2b', 'Uniform random gene sets (10,000), empirical one-sided p',
  'Mean number of TFs per target (and % of targets with >= 3 TFs) of the genes regulated by the five fitness-enhancing TFs '
  '(RbsR, GcvA, Fis, SlyA, RcsB) vs. 10,000 random gene sets of the same size from the RegulonDB 11.2 network; each '
  'fitness-enhancing TF and each global TF removed in turn. No random set reached the observed values (p < 0.0001).',
  _dlt_table(_starts[0]), 'input_DLT/tablas_supp.xlsx (D. Ledezma-Tejeida), Table S1'),
 ('S20w_TF_redundancy_networks', 'R2-M2b', 'Degree-preserving network randomization (curveball, 10,000 networks), empirical one-sided p',
  'Same statistics vs. 10,000 randomizations of the network that keep the number of targets of each TF and of regulators '
  'of each gene; leave-one-out of the fitness-enhancing TFs and removal of the seven global TFs.',
  _dlt_table(_starts[1]), 'input_DLT/tablas_supp.xlsx (D. Ledezma-Tejeida), Table S2'),
]

out = os.path.join(os.path.dirname(HERE), 'Table_S20_statistical_tests.xlsx')
with pd.ExcelWriter(out) as w:
    idx = pd.DataFrame([{'Sheet': s, 'Reviewer comment': rc, 'Test': te, 'Description': de, 'Rows': len(t), 'Source file (R1_analyses/output)': src}
                        for s, rc, te, de, t, src in S])
    pd.DataFrame([['Supplementary Table S20. Statistical tests added in the revision. Each sheet lists the reviewer comment it '
                   'addresses, the test and the source file of the analysis scripts (github.com/utrillalab/Fitness_analysis).']]
                 ).to_excel(w, sheet_name='Index', index=False, header=False)
    idx.to_excel(w, sheet_name='Index', index=False, startrow=2)
    for s, rc, te, de, t, src in S:
        t = t.copy()
        if isinstance(t.columns, pd.MultiIndex):
            t.columns = [' | '.join(str(x) for x in c if str(x).strip() and not str(x).startswith('Unnamed')) for c in t.columns]
        pd.DataFrame([[f'{s}. {de}'], [f'Reviewer comment: {rc}. Test: {te}. Source: {src}']]).to_excel(
            w, sheet_name=s[:31], index=False, header=False)
        t.to_excel(w, sheet_name=s[:31], index=False, startrow=3)
print(idx[['Sheet', 'Reviewer comment', 'Rows']].to_string(index=False))
print('core regulators:', x, 'of', n, f'fold {(x / n) / (K / N):.2f}', f'p {stats.hypergeom.sf(x - 1, N, K, n):.1e}')
print('saved', out)
