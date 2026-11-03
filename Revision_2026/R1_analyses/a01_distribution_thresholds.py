"""R1-M2 / R1-m2: fitness distributions with mode, mean and category cut-offs.
R1-M1: proteome fractions with genes lacking fitness data excluded vs treated as essential.
Sensitivity of proteome fractions to the mean-effect and essential cut-offs (also R2-M2a)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from common import *

fit, names = load_fitness()
mass = load_mass()
df = long_table()

# ---- 1. distribution summary per condition ----------------------------------------------------
LO, HI = -4, 2                                     # plotting window; values outside are pooled in edge bars
BINS = np.arange(LO, HI + 0.001, 0.05)


def summarize(f):
    kde = gaussian_kde(f, bw_method=0.05)
    grid = np.linspace(-3, 2, 2001)
    mode = grid[np.argmax(kde(grid))]
    mode = 0.0 if abs(mode) < 0.005 else mode
    core = f[np.abs(f - mode) < 0.5]
    sd = 1.4826 * np.median(np.abs(core - np.median(core)))
    return mode, sd


def panel(ax, f, title, color):
    mode, sd = summarize(f)
    ax.axvspan(-0.5, 0.5, color='#d9d9d9', alpha=0.5, lw=0, zorder=0)
    inside = f[(f > LO) & (f < HI)]
    ax.hist(inside, bins=BINS, color=color, zorder=2)
    w = BINS[1] - BINS[0]
    for edge, n in [(LO, (f <= LO).sum()), (HI - w, (f >= HI).sum())]:   # edge bars, hatched
        if n:
            ax.bar(edge, n, width=w, align='edge', color='white', edgecolor=color, hatch='////', lw=0.6, zorder=2)
    ax.axvline(-2, color='k', ls=':', lw=0.9, zorder=3)
    for x in (-0.5, 0.5):
        ax.axvline(x, color='k', ls='--', lw=0.8, zorder=3)
    ax.axvline(mode, color='#2166ac', lw=1.4, zorder=4)
    ax.axvline(f.mean(), color='#b2182b', lw=1.4, zorder=4)
    pct_me = 100 * np.mean((f >= -0.5) & (f <= 0.5))
    ax.text(0.03, 0.96, f'mode {mode:.2f}   mean {f.mean():+.2f}\nσ peak {sd:.2f}   mean-effect {pct_me:.0f}%',
            transform=ax.transAxes, ha='left', va='top', fontsize=8.5,
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.85), zorder=5)
    ax.set_yscale('log'); ax.set_title(title, fontsize=11)
    ax.set_xlim(LO - 0.05, HI + 0.05)
    ax.set_xticks([-4, -3, -2, -1, -0.5, 0, 0.5, 1, 2])
    ax.set_xticklabels(['≤−4', '−3', '−2', '−1', '', '0', '', '1', '≥2'])
    ax.spines[['top', 'right']].set_visible(False)
    return mode, sd


rows = []
fig, axes = plt.subplots(3, 4, figsize=(16, 10.5), sharex=True)
for ax, c in zip(axes.flat, CONDS):
    f = fit[c].dropna().values
    mode, sd = panel(ax, f, c, '#8c8c8c')
    cats = pd.Series(classify(f)).value_counts(normalize=True) * 100
    rows.append({'condition': c, 'n_genes': len(f), 'mode': round(mode, 3), 'median': round(np.median(f), 3),
                 'mean': round(f.mean(), 3), 'skewness': round(pd.Series(f).skew(), 2),
                 'robust_sd_central_peak': round(sd, 3),
                 **{f'%{k}': round(cats.get(k, 0), 1) for k in CATS}})
pooled = fit.values.ravel(); pooled = pooled[~np.isnan(pooled)]
panel(axes.flat[-1], pooled, 'All 11 conditions pooled', '#3a3a3a')
fig.supxlabel('Fitness (log2 ratio of final/initial mutant abundance)', fontsize=12, y=0.045)
fig.supylabel('Number of genes (log scale)', fontsize=12)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
handles = [Line2D([], [], color='#2166ac', lw=1.6, label='Mode'),
           Line2D([], [], color='#b2182b', lw=1.6, label='Mean'),
           Patch(fc='#d9d9d9', ec='k', ls='--', label='Mean-effect range (−0.5 to 0.5)'),
           Line2D([], [], color='k', ls=':', lw=1, label='Essential cut-off (−2)'),
           Patch(fc='white', ec='#3a3a3a', hatch='////', label='Values beyond the axis (≤ −4 or ≥ 2)')]
fig.legend(handles=handles, loc='lower center', ncol=5, frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0))
fig.tight_layout(rect=(0.01, 0.07, 1, 1))
fig.savefig(f'{OUT}/a01_fitness_distributions.png', dpi=200)
fig.savefig(f'{OUT}/a01_fitness_distributions.svg')
dist = pd.DataFrame(rows).set_index('condition')
dist.to_csv(f'{OUT}/a01_distribution_summary.csv')
print(dist.to_string())


# ---- 2. sensitivity of proteome fractions to cut-offs ----------------------------------------
def fractions(d, lo, hi, ess):
    d = d.assign(cat=classify(d.fitness, lo, hi, ess))
    tot = d.groupby('condition').mass_fg.sum()
    return (d.groupby(['condition', 'cat']).mass_fg.sum().unstack().div(tot, axis=0) * 100).reindex(columns=CATS)


sens = []
for lo, hi in [(-0.3, 0.3), (-0.5, 0.5), (-0.7, 0.7), (-1.0, 1.0)]:
    for ess in [-1.5, -2.0, -3.0]:
        fr = fractions(df, lo, hi, ess)
        for c in CONDS:
            sens.append({'mean_effect_bounds': f'+-{hi}', 'essential_cutoff': ess, 'condition': c,
                         **{k: round(fr.loc[c, k], 1) if pd.notna(fr.loc[c, k]) else 0 for k in CATS}})
sens = pd.DataFrame(sens)
sens.to_csv(f'{OUT}/a01_threshold_sensitivity.csv', index=False)
summ = sens.groupby(['mean_effect_bounds', 'essential_cutoff'])[CATS].agg(['min', 'median', 'max']).round(1)
summ.to_csv(f'{OUT}/a01_threshold_sensitivity_summary.csv')
print('\nProteome fraction (%) across 11 conditions under alternative cut-offs (min / median / max):')
print(summ.to_string())

# ---- 3. genes without fitness data: excluded (as in Fig. 1C) vs treated as essential (R1-M1) --
missing = sorted(b for b in set(mass.index) - set(fit.index) if isinstance(b, str))
extra = pd.concat([pd.DataFrame({'bnum': missing, 'condition': c, 'fitness': -10.0,
                                 'mass_fg': mass.loc[missing, c].values}) for c in CONDS])
excl = fractions(df, -0.5, 0.5, -2)
incl = fractions(pd.concat([df, extra]), -0.5, 0.5, -2)
m1 = pd.concat({'excluded (Fig. 1C)': excl, 'treated as essential': incl}, axis=1).round(1)
m1.to_csv(f'{OUT}/a01_missing_genes_scenarios.csv')
print(f'\n{len(missing)} proteomics genes lack fitness data. Their median share of the proteome: '
      f'{(mass.loc[missing].sum() / mass.sum() * 100).median():.1f}%')
print(m1.to_string())
