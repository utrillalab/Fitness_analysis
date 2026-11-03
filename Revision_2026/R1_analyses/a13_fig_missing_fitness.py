"""R1-M1 / R2-M1: supplementary panel showing the proteome share of proteins without fitness data.

Denominator = all proteins measured by Schmidt et al. (2016) in each condition, i.e. the 1,952 genes with
fitness data (Fig. 1C) plus the measured proteins absent from the RB-TnSeq dataset (no insertions or too few
reads; 243 of 395 (393 after the gatR/wbbL fix in common.py) have no Keio deletion mutant). Mean-effect fractions are reported for both denominators.
Colours: colour-blind-safe palette used for Fig. S4."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *

PAL = {'essential': '#b2182b', 'important': '#e08214', 'mean-effect': '#969696', 'enhancing': '#2166ac',
       'no fitness data': '#b2182b'}
ORDER = ['LB', 'Glycerol', 'Galactose', 'Mannose', 'Pyruvate', 'Glucosamine', 'Fructose', 'Glucose',
         'Xylose', 'Succinate', 'Acetate']                     # Fig. 1C order (glucosamine HCl after pyruvate)

fit, _ = load_fitness()
mass = load_mass()
df = long_table()
missing = sorted(b for b in set(mass.index) - set(fit.index) if isinstance(b, str))

rows = []
for c in ORDER:
    d = df[df.condition == c]
    m_missing = mass.loc[missing, c].sum()
    tot_all = d.mass_fg.sum() + m_missing
    r = {'condition': c, 'n_missing_with_mass>0': int((mass.loc[missing, c] > 0).sum())}
    for k in CATS:
        r[k] = 100 * d.loc[d.category == k, 'mass_fg'].sum() / tot_all
    r['no fitness data'] = 100 * m_missing / tot_all
    r['mean-effect_%_Fig1C_denominator'] = 100 * d.loc[d.category == 'mean-effect', 'mass_fg'].sum() / d.mass_fg.sum()
    rows.append(r)
res = pd.DataFrame(rows).set_index('condition')
res.round(2).to_csv(f'{OUT}/a13_missing_fitness_fractions.csv')
print(res.round(1).to_string())
print('\nmean-effect, Fig. 1C denominator: %.1f-%.1f%%; all measured proteins: %.1f-%.1f%%; '
      'no-fitness-data share: %.1f-%.1f%% (median %.1f%%)' % (
          res['mean-effect_%_Fig1C_denominator'].min(), res['mean-effect_%_Fig1C_denominator'].max(),
          res['mean-effect'].min(), res['mean-effect'].max(),
          res['no fitness data'].min(), res['no fitness data'].max(), res['no fitness data'].median()))

fig, ax = plt.subplots(figsize=(8.5, 4.6))
x = np.arange(len(ORDER)); bottom = np.zeros(len(ORDER))
stack = ['mean-effect', 'enhancing', 'important', 'essential', 'no fitness data']
for k in stack:
    v = res.loc[ORDER, k].values
    kw = dict(color=PAL[k], edgecolor='white', lw=0.6)
    if k == 'no fitness data':
        kw = dict(color='white', edgecolor=PAL[k], hatch='////', lw=0.8)
    ax.bar(x, v, 0.75, bottom=bottom, label=k, **kw)
    bottom += v
for i, c in enumerate(ORDER):
    ax.text(i, res.at[c, 'mean-effect'] / 2, f"{res.at[c, 'mean-effect']:.0f}", ha='center', va='center',
            fontsize=8, color='black')
ax.set_xticks(x); ax.set_xticklabels(ORDER, rotation=45, ha='right')
ax.set_ylabel('Measured proteome (%)'); ax.set_ylim(0, 100)
ax.spines[['top', 'right']].set_visible(False)
h, l = ax.get_legend_handles_labels()
ax.legend(h[::-1], l[::-1],
          frameon=False, fontsize=8, loc='upper left', bbox_to_anchor=(1.0, 1.0))
ax.set_title('Proteome mass by fitness category, including proteins without fitness data', fontsize=10)
fig.tight_layout()
for ext in ('png', 'svg', 'pdf'):
    fig.savefig(f'{OUT}/a13_missing_fitness.{ext}', dpi=300 if ext == 'png' else None)
