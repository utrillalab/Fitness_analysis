"""Compare the ViSEAGO GO-enrichment results of a15_viseago.R.

  repo vs original     : reproducibility (same inputs; differences = newer GO annotation / package versions)
  original vs corrected: effect of correcting the fitness values (Table S1 misalignment, glucosamine HCl)
  corrected vs no287   : R1-M1 sensitivity, essential category with and without the 287 essential genes without
                         fitness data
A term counts as enriched in a category when its p-value < 0.01 (the cut-off of the original analysis).
Output: output/a15b_viseago_comparison.csv, output/a15b_essential_terms_with_vs_without_287.csv
"""
import os
import pandas as pd
from common import *

CONDS = ['Glucose', 'LB', 'Glucosamine', 'Glycerol', 'Xylose', 'Mannose', 'Galactose', 'Fructose', 'Acetate', 'Pyruvate',
         'Succinate']
CATS = ['Essential', 'Important', 'No_effect', 'Fitness_enhancing']
SRC = {'repo': os.path.join(os.path.dirname(HERE), 'original_submission', 'Viseago_Res'),  # results as first submitted
       'original': os.path.join(OUT, 'viseago', 'original'),
       'corrected': os.path.join(OUT, 'viseago', 'corrected'),
       'no287': os.path.join(OUT, 'viseago', 'no287')}


def load(src, cond):
    d = pd.read_csv(os.path.join(SRC[src], f'BP_Results_{cond}.xls'), sep='\t')
    out = {}
    for c in CATS:
        p = pd.to_numeric(d.get(f'{c}.pvalue'), errors='coerce')
        out[c] = d.loc[p < 0.01, ['GO.ID', 'term']].assign(p=p[p < 0.01]).set_index('GO.ID')
    return out


def jac(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 1.0


rows, ess = [], []
for cond in CONDS:
    r = {s: load(s, cond) for s in SRC}
    for c in CATS:
        for a, b in [('repo', 'original'), ('original', 'corrected'), ('corrected', 'no287')]:
            A, B = r[a][c], r[b][c]
            rows.append({'condition': cond, 'category': c, 'comparison': f'{a} vs {b}', f'n_{a}': len(A), f'n_{b}': len(B),
                         'n_a': len(A), 'n_b': len(B), 'shared': len(set(A.index) & set(B.index)),
                         'jaccard': round(jac(A.index, B.index), 3),
                         'top10_kept': len(set(A.sort_values('p').index[:10]) & set(B.index))})
    W, N = r['corrected']['Essential'], r['no287']['Essential']
    for go in sorted(set(W.index) | set(N.index)):
        ess.append({'condition': cond, 'GO.ID': go, 'term': (W.term.get(go) if go in W.index else N.term.get(go)),
                    'p_with_287': W.p.get(go), 'p_without_287': N.p.get(go)})
res = pd.DataFrame(rows)[['condition', 'category', 'comparison', 'n_a', 'n_b', 'shared', 'jaccard', 'top10_kept']]
res.to_csv(os.path.join(OUT, 'a15b_viseago_comparison.csv'), index=False)
pd.DataFrame(ess).to_csv(os.path.join(OUT, 'a15b_essential_terms_with_vs_without_287.csv'), index=False)
summ = res.groupby(['comparison', 'category']).agg(n_a=('n_a', 'sum'), n_b=('n_b', 'sum'), shared=('shared', 'sum'),
                                                  jaccard_median=('jaccard', 'median'), top10_kept=('top10_kept', 'sum'))
print(summ.to_string())
