"""Figure of the fitness ~ log10(protein mass) regressions (R1-M4a, M4d; Table S20g).

Each panel: all gene x condition points of one group (11 conditions pooled) and the fitted line from a linear
regression with condition as covariate (slope = within-condition change in fitness per 10-fold increase in
mass; gene-clustered SE as in a02_mass_by_category.py). The line is drawn through the group means."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from common import *

PAL = {'all genes': '#4d4d4d', 'essential': '#b2182b', 'important': '#e08214', 'mean-effect': '#969696',
       'enhancing': '#2166ac'}

df = long_table()
df = df[df.mass_fg > 0].copy()
df['logm'] = np.log10(df.mass_fg)


def ols_cluster(y, X, groups):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ beta
    XtX_inv = np.linalg.pinv(X.T @ X)
    meat = np.zeros((X.shape[1], X.shape[1]))
    for g in np.unique(groups):
        idx = groups == g
        s = X[idx].T @ e[idx]
        meat += np.outer(s, s)
    G = len(np.unique(groups)); n, k = X.shape
    V = XtX_inv @ meat @ XtX_inv * G / (G - 1) * (n - 1) / (n - k)
    return beta, np.sqrt(np.diag(V)), G


groups = [('all genes', df)] + [(k, df[df.category == k]) for k in CATS]
fig, axes = plt.subplots(1, 5, figsize=(16, 3.8), sharex=True, sharey=True)
rows = []
for ax, (name, d) in zip(axes, groups):
    cond = pd.get_dummies(d.condition, drop_first=True, dtype=float)
    X = np.column_stack([np.ones(len(d)), d.logm.values, cond.values])
    beta, se, G = ols_cluster(d.fitness.values, X, d.bnum.values)
    slope, s = beta[1], se[1]
    p = 2 * stats.t.sf(abs(slope / s), G - 1)
    rows.append({'group': name, 'n_points': len(d), 'n_genes': G, 'slope': slope, 'SE': s, 'p': p})
    ax.scatter(d.logm, d.fitness.clip(-8, 8), s=3, alpha=0.15, color=PAL[name], lw=0, rasterized=True)
    xs = np.linspace(d.logm.quantile(0.01), d.logm.quantile(0.99), 50)
    ys = d.fitness.mean() + slope * (xs - d.logm.mean())
    ax.plot(xs, ys, color='black', lw=2)
    ptxt = 'p < 10⁻⁴' if p < 1e-4 else f'p = {p:.1g}'
    ax.set_title(f'{name} (n = {G} genes)', fontsize=10)
    ax.text(0.03, 0.04, f'slope {slope:+.3f} per 10× mass\n{ptxt}', transform=ax.transAxes, fontsize=8.5,
            va='bottom', bbox=dict(fc='white', ec='none', alpha=0.8))
    ax.axhline(0, color='k', lw=0.5, ls=':')
    ax.spines[['top', 'right']].set_visible(False)
axes[0].set_ylabel('Fitness'); axes[0].set_xlim(-7, 1.5)
fig.supxlabel('log10 protein mass (fg per cell)', fontsize=10)
fig.tight_layout()
for ext in ('png', 'svg', 'pdf'):
    fig.savefig(f'{OUT}/Fig_S1_regression.{ext}', dpi=300 if ext == 'png' else None)
res = pd.DataFrame(rows)
res.to_csv(f'{OUT}/a19_regression_slopes.csv', index=False)
print(res.round(4).to_string(index=False))
