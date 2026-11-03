"""Revised Supplementary Table S14 (R1-m1: transport sentence of the Abstract): GO terms enriched among the core
mean-effect genes (mean-effect in all 46 conditions), from a29_table_S14.R, with the proteome fraction of each term in
glucose computed as in Fig05_p1_ViSEAGO_listag_noeff.ipynb (summed glucose mass of the term genes / total glucose mass of
Infiles/allgen4691_fitfg_newname.csv, %) and the comparison with the original table.
Output: ../Table_S14_revised.xlsx, output/a29_table_S14_old_vs_new.csv
"""
import os
import pandas as pd
from common import *

ROOT = os.path.dirname(HERE)
new = pd.read_csv(f'{OUT}/viseago/core_mean_effect/BP_Results_core_mean_effect.xls', sep='\t')
old = pd.read_excel(os.path.join(ROOT, 'external', 'Supp_Tables_1-19_Hidalgo_et_al.xlsx'), sheet_name='Table S14', header=1)
allg = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'allgen4691_fitfg_newname.csv'))
tot = allg['Glucose'].sum()
core = pd.read_excel(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'core_no_effect.xlsx'))
cb = set(core.bnumber)

COND = ['Glucose', 'LB', 'Glucosamine', 'Glycerol', 'Xylose', 'Mannose', 'Galactose', 'Fructose', 'Acetate', 'Pyruvate',
        'Succinate']
core_frac = (100 * allg[allg.bnum.isin(cb)][COND].sum() / allg[COND].sum()).round(1)


def mass_pct(symbols):
    genes = [s.strip() for s in str(symbols).split(';') if s.strip()]
    return 100 * allg.loc[allg.gene.isin(genes), 'Glucose'].sum() / tot


p = 'Core_mean_effect.'
new['Proteome fraction (glucose, %)'] = new[p + 'Significant_genes_symbol'].map(mass_pct).round(3)
new = new.sort_values('Proteome fraction (glucose, %)', ascending=False).reset_index(drop=True)
new['In original Table S14'] = new['GO.ID'].isin(old['GO.ID']).map({True: 'yes', False: 'new'})
lost = old[~old['GO.ID'].isin(new['GO.ID'])][['GO.ID', 'term']]

cmp = pd.concat([new[['GO.ID', 'term', p + 'pvalue', 'Proteome fraction (glucose, %)', 'In original Table S14']],
                 lost.assign(**{'In original Table S14': 'only in original'})], ignore_index=True)
cmp.to_csv(f'{OUT}/a29_table_S14_old_vs_new.csv', index=False)

note = ('Supplementary Table S14 (revised). GO biological processes enriched among the 1,946 core mean-effect genes '
        '(mean-effect in all 46 conditions; topGO elim + Fisher, p < 0.01, nodeSize = 20; background: all E. coli K-12 '
        'genes in org.EcK12.eg.db), sorted by their proteome fraction in glucose (summed protein mass of the term genes, '
        '% of the measured proteome; Schmidt et al. 2016).')
with pd.ExcelWriter(os.path.join(ROOT, 'Table_S14_revised.xlsx')) as w:
    pd.DataFrame([[note]]).to_excel(w, sheet_name='Table S14', index=False, header=False)
    new.to_excel(w, sheet_name='Table S14', index=False, startrow=2)
    core_frac.rename('Core mean-effect proteome fraction (%)').to_frame().to_excel(w, sheet_name='Core fraction by condition')
print(f'terms: {len(new)} (original {len(old)}); new: {(new["In original Table S14"] == "new").sum()}; '
      f'lost: {len(lost)}')
print(new[['GO.ID', 'term', p + 'genes_frequency', p + 'pvalue', 'Proteome fraction (glucose, %)', 'In original Table S14']]
      .to_string(index=False))
print('lost:', lost.to_string(index=False))
print('core mean-effect fraction by condition (%):', core_frac.to_dict())
