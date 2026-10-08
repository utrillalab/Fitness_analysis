"""Shared data loading for the LTEE analysis (Supplementary Figure 4).

Fitness: mean of the two replicates of each condition in Input_Data/fit_organism_Keio.tsv (Price et al. 2018);
glucose = set1IT003/IT004. Categories as in the manuscript: essential <= -2; important -2 to -0.5;
mean-effect -0.5 to 0.5; enhancing > 0.5.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
INPUT = os.path.join(HERE, 'Input_Data')
OUT = os.path.join(HERE, 'Output')
os.makedirs(OUT, exist_ok=True)

CONDS = ['LB', 'Glycerol', 'Galactose', 'Mannose', 'Glucosamine', 'Pyruvate',
         'Fructose', 'Glucose', 'Xylose', 'Succinate', 'Acetate']
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
