"""Fig. 2 with the correct glucosamine condition (D-glucosamine hydrochloride, IT063/IT064).

Reproduces cell 31 of Code_Notebooks/Fig01_Panels_and_Boxplots.ipynb (fitness vs protein mass contribution of every
gene with fitness and proteomics data in the 11 shared conditions; same colours, cut-offs, axis limits and panel
order), with two changes: the glucosamine panel uses D-glucosamine hydrochloride (column 26 of the fitness averages;
the original used column 25, N-acetylglucosamine) and the panels are labelled A-K, as in the manuscript (the notebook
labels them C-M). Data loading as in fig1cd.py.
Outputs (R1_analyses/output/): Fig02_glucosamineHCl.{svg,pdf,png} (npj style, npj_style.py, 28 Oct 2026); fig2_glucosamine_NAG_vs_HCl.csv (genes per category)
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('FA_REPO', os.path.dirname(HERE))
OUT = os.path.join(HERE, 'R1_analyses', 'output')

table = pd.read_excel(os.path.join(REPO, 'Input_Data', '41587_2016_BFnbt3418_MOESM18_ESM_Mod.xlsx'), sheet_name='Table S6', skiprows=[0, 1])
pc = table.iloc[:, [0, 2, 5] + list(range(7, 29)) + [74, 75, 76, 77]]
pm = table.iloc[:, [0, 2, 5] + list(range(29, 51)) + [74, 75, 76, 77]]
pm.columns = pc.columns
pm = pm.drop(labels=2123, axis=0).drop_duplicates('Bnumber').set_index('Bnumber')

ft = pd.read_csv(os.path.join(REPO, 'Input_Data', 'fit_organism_Keio.tsv'), sep='\t')
fa = pd.DataFrame(index=ft.sysName.values)
fa['bnumber'] = ft.sysName.values
for c in range(5, 95, 2):
    fa[ft.columns[c]] = ft.iloc[:, [c, c + 1]].mean(axis=1).values
fa['LB'] = ft.iloc[:, [112, 131]].mean(axis=1).values
NAMES = ['bnumber', 'Glucose', 'Acetate', 'Glucosamine', 'Glycerol', 'Pyruvate', 'Xylose', 'Mannose', 'Galactose',
         'Succinate', 'Fructose', 'LB']
assert fa.columns[25].startswith('set1IT061 NAG') and fa.columns[26].startswith('set1IT063 D-Glucosamine Hydrochloride')
SEL = {'NAG': [0, 1, 12, 25, 19, 15, 4, 27, 5, 16, 2, 46], 'HCl': [0, 1, 12, 26, 19, 15, 4, 27, 5, 16, 2, 46]}


def selected(key):
    sf = fa.iloc[:, SEL[key]].copy(); sf.columns = NAMES
    return sf


def category(f):
    return np.select([f > 0.5, f < -2, (f > -2) & (f < -0.5)], ['enhancing', 'essential', 'important'], 'mean-effect')


# npj style (28 Oct 2026): Arial 7 pt at the printed width (180 mm), category colours of npj_style.py (as in Fig. 4),
# no title inside the figure (it repeats the legend title), bold lower-case letters
from npj_style import *
COL = {k: (PAL[k], 0.10 if k == 'mean-effect' else 0.55) for k in PAL}
sf = selected('HCl')
shared = sorted(set(sf.index) & set(pm.index))
fig = plt.figure(figsize=(DOUBLE, 135 * MM))
LEG = [('mean-effect', 'Mean-effect'), ('important', 'Important'), ('essential', 'Essential'), ('enhancing', 'Enhancing')]
fig.legend([plt.Line2D([], [], color=COL[k][0], marker='o', linestyle='none', markersize=4) for k, _ in LEG],
           [lab for _, lab in LEG], loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0))
axs = []
for n, cond in enumerate(NAMES[1:]):
    ax = plt.subplot(3, 4, n + 1); axs.append(ax)
    ax.set_title(cond)
    m = pd.to_numeric(pm.loc[shared, cond], errors='coerce'); f = sf.loc[shared, cond]
    ok = m > 0
    cat = pd.Series(category(f[ok].values), index=f[ok].index)
    for k in ['mean-effect', 'important', 'essential', 'enhancing']:     # draw grey first, as the many grey points are faint
        s_ = cat == k
        ax.plot(m[ok][s_], f[ok][s_], color=COL[k][0], alpha=COL[k][1], linestyle='none', marker='o', markersize=1.6,
                markeredgewidth=0, rasterized=True)
    ax.set_xlim(0.0000001, 33); ax.set_ylim(-8, 7); ax.set_xscale('log')
    ax.set_xticks([1e-7, 1e-5, 1e-3, 1e-1, 1e1])
    ax.set_xlabel('Protein mass (fg)'); ax.set_ylabel('Fitness')
plt.tight_layout(rect=(0, 0, 1, 0.95), h_pad=1.2, w_pad=1.5)
for n, ax in enumerate(axs):
    letter(fig, ax, chr(97 + n), dx=-0.005)
save(fig, os.path.join(OUT, 'Fig02_glucosamineHCl'))

rows = []
for key in ('NAG', 'HCl'):
    s = selected(key); m = pd.to_numeric(pm.loc[shared, 'Glucosamine'], errors='coerce'); f = s.loc[shared, 'Glucosamine']
    ok = m > 0
    rows.append(pd.Series(category(f[ok].values)).value_counts().rename(key))
cmp = pd.concat(rows, axis=1)
cmp.to_csv(os.path.join(OUT, 'fig2_glucosamine_NAG_vs_HCl.csv'))
print('Glucosamine panel, genes per category (NAG vs HCl):'); print(cmp.to_string())
