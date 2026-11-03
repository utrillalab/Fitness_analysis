"""Shared data loading for the Reviewer 1 analyses.

Fitness: mean of the two replicates of the SAME condition (as `fitness_averages` in
Final_Table_and_Proteome_Fractions.ipynb), with glucosamine taken from D-Glucosamine HCl
(IT063/IT064), the condition used by Schmidt et al. 2016 (the notebook used NAG, IT061/IT062).
Proteomics: Schmidt et al. 2016, Table S6, protein mass (fg/cell).
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('FA_REPO', os.path.dirname(os.path.dirname(HERE)))  # repository root
OUT = os.path.join(HERE, 'output')
os.makedirs(OUT, exist_ok=True)

CONDS = ['LB', 'Glycerol', 'Galactose', 'Mannose', 'Glucosamine', 'Pyruvate',
         'Fructose', 'Glucose', 'Xylose', 'Succinate', 'Acetate']
ME_CONDS = [c for c in CONDS if c != 'LB']
CATS = ['essential', 'important', 'mean-effect', 'enhancing']
COLORS = {'essential': '#b2182b', 'important': '#ef8a62', 'mean-effect': '#636363', 'enhancing': '#1b7837'}
FIT_COL = {'Glucose': 'set1IT003 D-Glucose (C)', 'Acetate': 'set1IT027 acetate (C)',
           'Glucosamine': 'set1IT063 D-Glucosamine Hydrochloride (C)', 'Glycerol': 'set1IT047 Glycerol (C)',
           'Pyruvate': 'set1IT037 pyruvate (C)', 'Xylose': 'set1IT011 D-Xylose (C)',
           'Mannose': 'set1IT065 D-Mannose (C)', 'Galactose': 'set1IT013 D-Galactose (C)',
           'Succinate': 'set1IT039 succinate (C)', 'Fructose': 'set1IT005 D-Fructose (C)', 'LB': 'LB'}


def classify(f, lo=-0.5, hi=0.5, ess=-2.0):
    """Manuscript rule; thresholds can be changed for sensitivity analyses."""
    f = np.asarray(f, dtype=float)
    out = np.where(f <= ess, 'essential',
          np.where(f < lo, 'important',
          np.where(f <= hi, 'mean-effect', 'enhancing')))
    return np.where(np.isnan(f), 'NA', out)


def load_fitness_all():
    """All 46 conditions (45 replicate pairs + LB), genes x conditions, plus gene names."""
    ft = pd.read_csv(os.path.join(REPO, 'Input_Data', 'fit_organism_Keio.tsv'), sep='\t')
    fa = pd.DataFrame(index=ft.sysName.values)
    for c in range(5, 95, 2):
        fa[ft.columns[c]] = ft.iloc[:, [c, c + 1]].mean(axis=1).values
    fa['LB'] = ft.iloc[:, [112, 131]].mean(axis=1).values
    names = pd.Series(ft.geneName.values, index=ft.sysName.values)
    return fa, names


def load_fitness():
    fa, names = load_fitness_all()
    return pd.DataFrame({c: fa[FIT_COL[c]] for c in CONDS}), names


def load_mass():
    t = pd.read_excel(os.path.join(REPO, 'Input_Data', '41587_2016_BFnbt3418_MOESM18_ESM_Mod.xlsx'),
                      sheet_name='Table S6', skiprows=[0, 1])
    pc = t.iloc[:, [0, 2, 5] + list(range(7, 29)) + [74, 75, 76, 77]]
    pm = t.iloc[:, [0, 2, 5] + list(range(29, 51)) + [74, 75, 76, 77]]
    pm.columns = pc.columns
    pm = pm.drop(labels=2123, axis=0).drop_duplicates('Bnumber').set_index('Bnumber')
    pm = pm[pm.index.notna()]
    # gatR and wbbL carry their old pseudogene b-numbers in Schmidt et al.; the fitness data use the current ones
    # (b4498, b4571). yqjK has no positive mass in any condition. With these fixes: 1,952 genes with fitness and
    # proteomics and 393 measured proteins without fitness, as in Table S1 (audit B3, 26 Oct 2026).
    pm = pm.rename(index={'b2090': 'b4498', 'b2031': 'b4571'})
    pm = pm[CONDS].apply(pd.to_numeric, errors='coerce')
    return pm[(pm > 0).any(axis=1)]


def load_flux():
    f = pd.read_csv(os.path.join(REPO, 'ME_model', 'ME_Output', 'Updated_ME_fluxes.csv')).set_index('genes')
    return f[ME_CONDS]


def long_table():
    """One row per gene x condition for the 1,952 genes with both fitness and proteomics."""
    fit, names = load_fitness()
    mass = load_mass()
    genes = sorted(set(fit.index) & set(mass.index))
    rows = []
    for c in CONDS:
        d = pd.DataFrame({'bnum': genes, 'gene': names.loc[genes].values, 'condition': c,
                          'fitness': fit.loc[genes, c].values, 'mass_fg': mass.loc[genes, c].values})
        rows.append(d)
    df = pd.concat(rows, ignore_index=True)
    df['category'] = classify(df.fitness)
    return df


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    m = len(p)
    running = 0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(1, running)
    return adj


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)[::-1]
    adj = np.empty_like(p)
    prev = 1
    for rank, i in enumerate(order):
        prev = min(prev, p[i] * m / (m - rank))
        adj[i] = prev
    return adj


def cliffs_delta(x, y):
    x, y = np.asarray(x), np.asarray(y)
    from scipy.stats import mannwhitneyu
    u = mannwhitneyu(x, y, alternative='two-sided').statistic
    return 2 * u / (len(x) * len(y)) - 1


def cohens_d(x, y):
    x, y = np.asarray(x), np.asarray(y)
    nx, ny = len(x), len(y)
    s = np.sqrt(((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / (nx + ny - 2))
    return (x.mean() - y.mean()) / s
