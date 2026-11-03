"""Fig. 1C and 1D with the correct glucosamine condition (D-glucosamine hydrochloride, IT063/IT064).

Reproduces cells 13, 17, 21, 29 and 30 of Code_Notebooks/Final_Table_and_Proteome_Fractions.ipynb, with one change:
in cell 13 the "Glucosamine" column was taken from position 25 of fitness_averages (set1IT061 NAG = N-acetylglucosamine)
instead of position 26 (set1IT063 D-Glucosamine Hydrochloride), the condition of the Schmidt et al. 2016 proteomics.
Everything else (genes, denominators, colours, order) is as in the notebook.

Outputs (R1_analyses/output/):
  Fig01_revised.{pdf,svg,png}: complete Fig. 1 (a-d) in the npj style (28 Oct 2026; npj_style.py)
  fig1cd_fractions_NAG_vs_HCl.csv  (old vs new values for every condition)
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('FA_REPO', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(HERE, 'R1_analyses', 'output')

# cell 12: proteomics (fg per cell)
table = pd.read_excel(os.path.join(REPO, 'Input_Data', '41587_2016_BFnbt3418_MOESM18_ESM_Mod.xlsx'), sheet_name='Table S6', skiprows=[0, 1])
pc = table.iloc[:, [0, 2, 5] + list(range(7, 29)) + [74, 75, 76, 77]]
pm = table.iloc[:, [0, 2, 5] + list(range(29, 51)) + [74, 75, 76, 77]]
pm.columns = pc.columns
pm0 = pm.drop(labels=2123, axis=0)
pm = pm0.drop_duplicates('Bnumber').set_index('Bnumber')

# cell 13: fitness averages (mean of the two replicates of each condition)
ft = pd.read_csv(os.path.join(REPO, 'Input_Data', 'fit_organism_Keio.tsv'), sep='\t')
fa = pd.DataFrame(index=ft.sysName.values)
fa['bnumber'] = ft.sysName.values
for c in range(5, 95, 2):
    fa[ft.columns[c]] = ft.iloc[:, [c, c + 1]].mean(axis=1).values
fa['LB'] = ft.iloc[:, [112, 131]].mean(axis=1).values
NAMES = ['bnumber', 'Glucose', 'Acetate', 'Glucosamine', 'Glycerol', 'Pyruvate', 'Xylose', 'Mannose', 'Galactose',
         'Succinate', 'Fructose', 'LB']
COLS = {'NAG (original notebook)': [0, 1, 12, 25, 19, 15, 4, 27, 5, 16, 2, 46],
        'Glucosamine HCl (corrected)': [0, 1, 12, 26, 19, 15, 4, 27, 5, 16, 2, 46]}
assert fa.columns[25].startswith('set1IT061 NAG') and fa.columns[26].startswith('set1IT063 D-Glucosamine Hydrochloride')

fluxes = pd.read_csv(os.path.join(REPO, 'ME_model', 'ME_Output', 'Updated_ME_fluxes.csv')).set_index('genes')
C1C = ['LB', 'Glycerol', 'Galactose', 'Mannose', 'Pyruvate', 'Glucosamine', 'Fructose', 'Glucose', 'Xylose', 'Succinate', 'Acetate']  # descending mean-effect (glucosamine HCl 53.1% < pyruvate 54.7%)
C1D = C1C[1:]


def classify(f):
    return np.select([f <= -2, (f >= -0.5) & (f <= 0.5), (f > -2) & (f < -0.5), f > 0.5],
                     ['essential', 'no_effect', 'important', 'fitness-enhancing'], 'huh')


def fractions(sel):
    sf = fa.iloc[:, sel].copy(); sf.columns = NAMES
    shared = sorted(set(pm.index) & set(sf.index))
    res = {}
    for c in C1C:
        # cell 17 / 29: genes with fitness and proteomics; fraction of their total mass per category
        fg = pd.to_numeric(pm.loc[shared, c], errors='coerce')
        cl = pd.Series(classify(sf.loc[shared, c].values), index=shared)
        tot = np.nansum(fg)
        r = {k: np.nansum(fg[cl == k]) / tot for k in ['no_effect', 'essential', 'important', 'fitness-enhancing']}
        # cell 21 / 30: mean-effect mass split by ME-model status (same genes, same denominator)
        if c in fluxes.columns:
            fl = fluxes[c].reindex(shared)
            ne = cl == 'no_effect'
            r['ME_absent'] = np.nansum(fg[ne & fl.isna()]) / tot
            r['ME_zero_flux'] = np.nansum(fg[ne & (fl == 0)]) / tot
            r['ME_nonzero_flux'] = np.nansum(fg[ne & fl.notna() & (fl != 0)]) / tot
        res[c] = r
    return pd.DataFrame(res).T


tabs = {k: fractions(v) for k, v in COLS.items()}
cmp = pd.concat(tabs, axis=1).round(4)
cmp.to_csv(os.path.join(OUT, 'fig1cd_fractions_NAG_vs_HCl.csv'))
print((cmp * 100).round(1).loc[['Glucosamine', 'Glucose', 'Acetate', 'LB']].T.to_string())
same = (cmp['NAG (original notebook)'].drop('Glucosamine') - cmp['Glucosamine HCl (corrected)'].drop('Glucosamine')).abs().max().max()
print('max difference in the other conditions:', same)
new = tabs['Glucosamine HCl (corrected)']


# ---- plots: complete Fig. 1 in the npj style (28 Oct 2026) ------------------------------------------------
# a, b redrawn from the data with the code of Fig01_Panels_and_Boxplots.ipynb (cells 27 and 29: glucose, all genes with
# fitness / all proteins, same bins) instead of the raster of the original figure; c, d as before. Arial 7 pt at the
# printed width (180 mm); category colours as in Figs. 2 and 4; mean-effect subsets of d in greys (colour-blind safe).
from npj_style import *


def plot_1a(ax):
    v = fa.iloc[:, 1]   # Glucose (set1IT001)
    ax.hist(v.dropna(), bins=[-8 + x / 4 for x in range(40)], color='#4a90c2', edgecolor='black', lw=0.4)
    ax.set_yscale('log'); ax.set_xlim(-6.5, 2)   # no gene below -5.5
    ax.set_xlabel('Fitness'); ax.set_ylabel('Gene counts')


def plot_1b(ax):
    v = pd.to_numeric(pm0['Glucose'], errors='coerce')   # protein_mass of cell 9
    ax.hist(v[v > 0], bins=[1e-7 * 10 ** (0.4 * x) for x in range(22)], color='#4a90c2', edgecolor='black', lw=0.4)
    ax.set_xscale('log'); ax.set_xlabel('Protein mass contribution (fg)'); ax.set_ylabel('Gene counts')


def plot_1c(ax):
    for k, c in enumerate(C1C):
        r = new.loc[c]; b = 0
        for key, cat in [('no_effect', 'mean-effect'), ('essential', 'essential'), ('important', 'important'),
                         ('fitness-enhancing', 'enhancing')]:
            ax.bar(k, r[key], bottom=b, color=PAL[cat], width=0.8, label=LABEL[cat] if k == 0 else '')
            b += r[key]
    ax.set_xticks(range(len(C1C))); ax.set_xticklabels(C1C, rotation=60, ha='right', rotation_mode='anchor')
    ax.set_ylabel('Proteome mass fraction'); ax.set_ylim(0, 1)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[::-1], l[::-1], loc='center left', bbox_to_anchor=(1.0, 0.5), handlelength=1, handleheight=1)


def plot_1d(ax):
    for k, c in enumerate(C1D):
        r = new.loc[c]
        ax.bar(k, r['ME_absent'], color=GREYS[0], label='Absent from\nME-model' if k == 0 else '')
        ax.bar(k, r['ME_zero_flux'], bottom=r['ME_absent'], color=GREYS[1], label='Zero flux' if k == 0 else '')
        ax.bar(k, r['ME_nonzero_flux'], bottom=r['ME_absent'] + r['ME_zero_flux'], color=GREYS[2],
               label='Non-zero flux' if k == 0 else '')
    ax.set_title('Mean-effect genes')
    ax.set_xticks(range(len(C1D))); ax.set_xticklabels(C1D, rotation=60, ha='right', rotation_mode='anchor')
    ax.set_ylabel('Proteome mass fraction')
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[::-1], l[::-1], loc='center left', bbox_to_anchor=(1.0, 0.5), handlelength=1, handleheight=1)


fig = plt.figure(figsize=(DOUBLE, 125 * MM))
gs = fig.add_gridspec(2, 2, left=0.08, right=0.86, bottom=0.14, top=0.96, hspace=0.55, wspace=0.55)
axs = [fig.add_subplot(gs[i, j]) for i in range(2) for j in range(2)]
for ax, fn in zip(axs, (plot_1a, plot_1b, plot_1c, plot_1d)):
    fn(ax)
for ax, s in zip(axs, 'abcd'):
    letter(fig, ax, s, dx=-0.01)
save(fig, os.path.join(OUT, 'Fig01_revised'))
print('saved Fig01_revised')
