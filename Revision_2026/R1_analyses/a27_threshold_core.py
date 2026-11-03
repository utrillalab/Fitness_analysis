"""R1-m10: sensitivity of the regulator enrichment to the recurrence threshold, for the candidates and for the core.

Candidates: fitness > 0.5 in >= k of the 46 conditions. Core: candidates that also have fitness > 0.5 and t > 4 in >= 1
of those conditions (option B, as a05). Regulators: RegulonDB TFs (NetworkRegulatorGene.tsv), hypergeometric test
against all genes with fitness data (as a03_m10 / a05 / a16). t: combined replicate t of Price et al. 2018 (as a05).
Thresholds k = 3-25; each set also without rpsG and rpsR (R2-m2).
Output: output/a27_threshold_core_sensitivity.csv (read by a16_table_S20.py, sheet S20p)
"""
import os
import numpy as np, pandas as pd
from scipy import stats
from common import *

P18 = os.path.join(os.path.dirname(HERE), 'Input_Data_Price2018')
tt = pd.read_csv(os.path.join(P18, 'fit_t.tab'), sep='\t').set_index('sysName')
ft = pd.read_csv(os.path.join(REPO, 'Input_Data', 'fit_organism_Keio.tsv'), sep='\t').set_index('sysName')
fa, names46 = load_fitness_all()
t46 = {}
for col in fa.columns:
    if col == 'LB':
        t46[col] = tt[['set2IT045 LB', 'set2IT069 LB']].sum(axis=1) / np.sqrt(3)
    else:
        j = list(ft.columns).index(col)
        pair = [ft.columns[j], ft.columns[j + 1]]
        t46[col] = tt[pair].sum(axis=1) / np.sqrt(3) if all(p in tt.columns for p in pair) else np.nan
t46 = pd.DataFrame(t46).reindex(fa.index)
n_enh = (fa > 0.5).sum(axis=1)
n_enh_t = ((fa > 0.5) & (t46 > 4)).sum(axis=1)

net = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'NetworkRegulatorGene.tsv'), sep='\t', comment='#', header=None)
lower = {n.lower(): n for n in names46.dropna().unique()}
regs = {lower[r.lower()] for r in set(net[1].dropna().astype(str)) if r.lower() in lower}
reg_b = set(names46.index[names46.isin(regs)])
N, K = len(fa), len(reg_b & set(fa.index))

# R2-m2: rpsG and rpsR, the only candidates essential in Keio, PEC and DEG (probable Tn-seq artefacts), removed from the
# sets (background unchanged, as in Table S20q)
art_b = set(names46.index[names46.isin(['rpsG', 'rpsR'])])
rows = []
for k in range(3, 26):
    cand = set(n_enh.index[n_enh >= k])
    core = {g for g in cand if n_enh_t[g] >= 1}
    for label, s in [('candidates', cand), ('core (t>4 in >=1)', core)]:
        for art, s2 in [('all genes', s), ('without rpsG and rpsR', s - art_b)]:
            x = len(s2 & reg_b)
            rows.append({'k': k, 'set': label, 'artefacts': art, 'genes': len(s2), 'regulators': x,
                         'fold': round((x / len(s2)) / (K / N), 2), 'p': stats.hypergeom.sf(x - 1, N, K, len(s2)),
                         'probable artefacts in set': ', '.join(sorted(names46[list(s & art_b)])) if art == 'all genes' else ''})
out = pd.DataFrame(rows)
out.to_csv(f'{OUT}/a27_threshold_core_sensitivity.csv', index=False, float_format='%.2g')
print(out[out.k.isin([5, 8, 10, 12, 15, 20])].to_string(index=False))
