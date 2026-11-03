"""Check that the 46-condition lists (108 fitness-enhancing genes, core mean-effect genes, Table S19)
match fitness values averaged with correctly paired replicates (notebook cell 13 / cell 15 logic)."""
import os, pandas as pd, numpy as np
repo = os.environ.get('FA_REPO', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ft = pd.read_csv(repo + '/Input_Data/fit_organism_Keio.tsv', sep='\t')
fa = pd.DataFrame(index=ft.sysName)
for c in range(5, 95, 2): fa[ft.columns[c]] = ft.iloc[:, [c, c+1]].mean(axis=1).values
fa['LB'] = ft.iloc[:, [112, 131]].mean(axis=1).values
enh = (fa > 0.5).sum(axis=1); noef = ((fa >= -0.5) & (fa <= 0.5)).sum(axis=1)
print('conditions:', fa.shape[1])
e108 = pd.read_excel(repo + '/Functional_Analysis/Infiles/108_enhancing_10Conditions.xlsx', skiprows=2, header=None).dropna(subset=[0])
e108.columns = ['bnumber', 'gene', 'product', 'n_cond']
calc = enh[enh >= 10]
print(f'fitness-enhancing in >=10 conditions: recomputed {len(calc)}, file {len(e108)}; '
      f'same genes: {set(calc.index) == set(e108.bnumber)}; same counts: {all(calc.loc[e108.bnumber].values == e108.n_cond.values)}')
core = pd.read_excel(repo + '/Functional_Analysis/Infiles/core_no_effect.xlsx', skiprows=1, header=None).dropna(subset=[1])
calc_core = noef[noef == 46]
print(f'core mean-effect (46/46): recomputed {len(calc_core)}, file {len(core)}; same genes: {set(calc_core.index) == set(core[1])}')
s19 = pd.read_excel(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'external', 'Supp_Tables_1-19_Hidalgo_et_al.xlsx'), sheet_name='Table S19')
print('Table S19 columns:', list(s19.columns)[:12])
s19 = s19.set_index('bnumber')
g = s19.index.intersection(fa.index)
ok = all((s19.loc[g, 'sum_of_mean-effect'] == noef.loc[g]).values) and all((s19.loc[g, 'sum_of_enhancing'] == enh.loc[g]).values)
print(f'Table S19 ({len(g)} genes) matches correct 46-condition counts: {ok}')
