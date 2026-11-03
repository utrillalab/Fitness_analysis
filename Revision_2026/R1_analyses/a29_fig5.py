"""Redrawn Fig. 5 (core mean-effect genes) with the revised Table S14 (a29_table_S14.R / .py).

  A  proteome mass fraction of the 1,946 core mean-effect genes (mean-effect in all 46 conditions) in each condition,
     relative to the mass of all genes with fitness data (as in Final_Table_and_Proteome_Fractions.ipynb, cell 31);
  B  the 12 enriched GO processes of the core with the largest proteome fraction in glucose (summed mass of the term
     genes / total measured proteome, %; as in Fig05_p1_ViSEAGO_listag_noeff.ipynb), with the number of core genes
     of each term; terms are nested (e.g. most transport terms are part of transmembrane transport), so their
     fractions are not additive.
Output: output/Fig05_revised.{svg,pdf,png}
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *

import sys; sys.path.insert(0, os.path.dirname(HERE)); from npj_style import *   # npj style (28 Oct 2026)
COL, INK = '#5f9ea0', '#222222'   # cadetblue, as in the original figure
ORDER = ['LB', 'Glycerol', 'Galactose', 'Mannose', 'Pyruvate', 'Glucosamine', 'Fructose', 'Glucose', 'Xylose',
         'Succinate', 'Acetate']

fa, _ = load_fitness_all()
allg = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'allgen4691_fitfg_newname.csv'))
core = set(pd.read_excel(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'core_no_effect.xlsx')).bnumber)
f = allg[allg.bnum.isin(fa.index)]
frac = f[f.bnum.isin(core)][ORDER].sum() / f[ORDER].sum()

p = 'Core_mean_effect.'
t = pd.read_excel(os.path.join(os.path.dirname(HERE), 'Table_S14_revised.xlsx'), header=2).head(12)
t['n'] = t[p + 'genes_frequency'].str.extract(r'\((\d+)/')[0].astype(int)

fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 76 * MM), gridspec_kw={'width_ratios': [1, 1.3], 'wspace': 1.15})
a.bar(range(len(ORDER)), frac.values, 0.7, color=COL)
a.set_xticks(range(len(ORDER)), ORDER, rotation=60, ha='right')
a.set_ylabel('Proteome mass fraction'); a.set_ylim(0, 0.16)
a.set_title('Core mean-effect genes\n(n = 1,946)')

y = np.arange(len(t))[::-1]
b.barh(y, t['Proteome fraction (glucose, %)'], 0.7, color=COL)
for yi, v, n in zip(y, t['Proteome fraction (glucose, %)'], t.n):
    b.text(v + 0.02, yi, f'{v:.2f} [{n}]', va='center', fontsize=6, color=INK)
b.set_yticks(y, t.term)
b.set_xlabel('Proteome mass fraction in glucose (%)'); b.set_xlim(0, 1.6)
b.set_title('Enriched GO processes of the core (glucose)')
fig.subplots_adjust(left=0.08, right=0.97, bottom=0.27, top=0.86)
letter(fig, a, 'a', dx=-0.01); letter(fig, b, 'b', dx=-0.01)
save(fig, f'{OUT}/Fig05_revised')
print(frac.round(3).to_dict())
print(t[['term', 'Proteome fraction (glucose, %)', 'n']].to_string(index=False))
