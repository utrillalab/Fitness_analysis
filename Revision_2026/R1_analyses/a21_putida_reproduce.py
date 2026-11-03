"""R1-M7a: independent reproduction of Supplementary Figure S5 (P. putida KT2440), from O. A. Aguilar-Vera's package
in ../external/Pp_FigureS5_share (notebook P.putida_Fitness_ExpressionTPM.ipynb, steps 1-3, 5 and 8).

Expression: GSE85475 (Bojanovic et al. 2017), TPM from raw counts (scripts/tpm_from_counts.py; rerun gives an identical
table), mean of control replicates C_1_0, C_2_0, C_3_0, log10(TPM + 1).
Fitness: Fitness Browser, experiment set6IT057 D-Glucose (C); manuscript categories.
Highly expressed: top 5% of log10(TPM + 1) among genes with fitness and expression.
COG test: top 5% vs remaining mean-effect genes (single-letter COGs), Fisher's exact test, BH within each level.
Checks every number against the package's results/ files.
Output: output/a21_putida_reproduction.csv
"""
import os
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from common import *

PP = os.path.join(os.path.dirname(HERE), 'external', 'Pp_FigureS5_share')
COG_GROUP = dict(zip('JAKLBDYVTMNZWUOXCGEFHIPQRS',
                     ['INFORMATION'] * 5 + ['CELLULAR PROCESSES'] * 11 + ['METABOLISM'] * 8 + ['POORLY CHARACTERIZED'] * 2))

tpm = pd.read_csv(os.path.join(PP, 'data', 'processed', 'GSE85475_TPM.tsv'), sep='\t')
tpm['logTPM'] = np.log10(tpm[['TPM_C_1_0', 'TPM_C_2_0', 'TPM_C_3_0']].mean(axis=1) + 1)
fit = pd.read_csv(os.path.join(PP, 'data', 'raw', 'fit_organism_Putida.tsv'), sep='\t')
d = fit[['locusId', 'set6IT057 D-Glucose (C)']].rename(columns={'set6IT057 D-Glucose (C)': 'fitness'}) \
    .merge(tpm[['locus_tag', 'logTPM']], left_on='locusId', right_on='locus_tag')
d['category'] = classify(d.fitness)
cut = d.logTPM.quantile(0.95)
d['top'] = d.logTPM >= cut
sel = d[(d.category == 'mean-effect') & d.top]

ft = pd.read_table(os.path.join(PP, 'data', 'raw', 'GCA_000007565_feature_table.txt'))
cogs = pd.read_table(os.path.join(PP, 'data', 'raw', 'COGs.csv')).rename(columns={'bnum': 'locus_tag'})
genes = ft.loc[ft['# feature'] == 'gene', ['locus_tag']].drop_duplicates().merge(cogs[['locus_tag', 'func']], how='left')
me = d[d.category == 'mean-effect'].merge(genes, how='left', on='locus_tag')
me1 = me[me.func.notna() & ~me.func.str.contains(',', na=False)].assign(group=lambda x: x.func.map(COG_GROUP))


def enrich(col):
    rows = []
    a_tot, b_tot = me1.top.sum(), (~me1.top).sum()
    for k in sorted(me1[col].dropna().unique()):
        a = ((me1[col] == k) & me1.top).sum(); b = ((me1[col] == k) & ~me1.top).sum()
        o, p = fisher_exact([[a, a_tot - a], [b, b_tot - b]])
        rows.append({'category': k, 'top_n': a, 'rest_n': b, 'odds_ratio': o, 'p': p})
    r = pd.DataFrame(rows); r['FDR'] = bh(r.p.values)
    return r


res = pd.concat([enrich('group').assign(level='COG group'), enrich('func').assign(level='COG letter')], ignore_index=True)
ref = pd.read_csv(os.path.join(PP, 'results', 'TableSX_COG_enrichment.tsv'), sep='\t')
m = res.merge(ref, on=['level', 'category'], suffixes=('', '_omar'))
tabC = pd.read_csv(os.path.join(PP, 'results', 'FigureS5_tableC.tsv'), sep='\t')

checks = [
    ('genes with fitness and expression', len(d), 4477),
    ('expression cutoff log10(TPM+1)', round(cut, 3), 2.400),
    ('genes above cutoff', int(d.top.sum()), 224),
    ('mean-effect genes above cutoff (Table C)', len(sel), len(tabC)),
    ('same Table C genes', int(set(sel.locusId) == set(tabC.locusId)), 1),
    ('mean-effect genes (all)', int((d.category == 'mean-effect').sum()), 4210),
    ('mean-effect genes with one COG letter', len(me1), 3315),
    ('COG tests matching (counts, p within 1e-6)', int(((m.top_n == m.top_n_omar) & (m.rest_n == m.rest_n_omar)
                                                       & np.isclose(m.p, m.p_omar, rtol=1e-3, atol=1e-6)).sum()), len(ref)),
]
out = pd.DataFrame(checks, columns=['quantity', 'reproduced', 'package'])
out['match'] = out.reproduced == out.package
out.to_csv(os.path.join(OUT, 'a21_putida_reproduction.csv'), index=False)
print(out.to_string(index=False))
print('\nCategories above the cutoff:', d[d.top].category.value_counts().to_dict())
print('All genes by category:', d.category.value_counts().to_dict(), f'(mean-effect {100 * (d.category == "mean-effect").mean():.1f}%)')
print('\nCOG tests with FDR < 0.05:\n', res[res.FDR < 0.05].round(4).to_string(index=False))
