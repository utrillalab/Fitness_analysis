"""Protein mass vs fitness category.
R1-m4 : observed mass fraction of each category vs a size-matched random-gene null.
R1-M4b: effect sizes (Cliff's delta, Cohen's d on log10 mass) for each category vs mean-effect.
R1-M4c: over/under-representation of categories across protein-mass quintiles.
R1-M4d: regression of fitness on log10 mass, per category, with condition fixed effects.
R1-m9 : fitness-enhancing vs mean-effect: location and spread of log10 mass."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from common import *

rng = np.random.default_rng(1)
df = long_table()
df = df[df.mass_fg > 0].copy()          # log scale; zero-abundance proteins dropped (reported below)
df['logm'] = np.log10(df.mass_fg)
print('gene x condition rows with mass > 0:', len(df))

# ---- m4: permutation null for mass fractions --------------------------------------------------
NPERM = 10000
perm = []
for c, d in df.groupby('condition'):
    m = d.mass_fg.values; tot = m.sum()
    for k in CATS:
        n = int((d.category == k).sum())
        if n == 0:
            continue
        obs = m[d.category.values == k].sum() / tot
        null = np.array([m[rng.choice(len(m), n, replace=False)].sum() for _ in range(NPERM)]) / tot
        p_hi = (1 + (null >= obs).sum()) / (NPERM + 1); p_lo = (1 + (null <= obs).sum()) / (NPERM + 1)
        perm.append({'condition': c, 'category': k, 'n_genes': n, 'gene_share_%': 100 * n / len(m),
                     'observed_mass_%': 100 * obs, 'null_mean_%': 100 * null.mean(),
                     'null_95CI_%': f'{100*np.percentile(null, 2.5):.1f}-{100*np.percentile(null, 97.5):.1f}',
                     'obs/expected': obs / null.mean(), 'p_two_sided': min(1, 2 * min(p_hi, p_lo))})
perm = pd.DataFrame(perm)
perm['p_BH'] = bh(perm.p_two_sided)
perm.round(4).to_csv(f'{OUT}/a02_m4_permutation_null.csv', index=False)
print('\nm4 - observed / expected mass fraction (size-matched null, 10,000 permutations):')
print(perm.pivot(index='condition', columns='category', values='obs/expected')[CATS].round(2).to_string())
print('significant (BH<0.05):', perm.assign(s=perm.p_BH < 0.05).groupby('category').s.sum().to_dict())

# ---- M4b: effect sizes vs mean-effect -----------------------------------------------------------
eff = []
for c, d in df.groupby('condition'):
    ref = d.loc[d.category == 'mean-effect', 'logm']
    for k in ['essential', 'important', 'enhancing']:
        x = d.loc[d.category == k, 'logm']
        if len(x) < 3:
            continue
        kw = stats.kruskal(x, ref)
        eff.append({'condition': c, 'category': k, 'n': len(x), 'n_mean_effect': len(ref),
                    'median_log10fg': x.median(), 'median_log10fg_mean_effect': ref.median(),
                    'fold_change_median_mass': 10 ** (x.median() - ref.median()),
                    'cliffs_delta': cliffs_delta(x, ref), 'cohens_d_log10': cohens_d(x, ref),
                    'kruskal_H': kw.statistic, 'p': kw.pvalue})
eff = pd.DataFrame(eff); eff['p_holm'] = holm(eff.p)
eff.round(4).to_csv(f'{OUT}/a02_M4b_effect_sizes.csv', index=False)
print('\nM4b - Cliff\'s delta (log10 mass, category vs mean-effect):')
print(eff.pivot(index='condition', columns='category', values='cliffs_delta').round(2).to_string())
print(eff.groupby('category')[['cliffs_delta', 'cohens_d_log10', 'fold_change_median_mass']].median().round(2).to_string())

# ---- M4c: enrichment across mass quintiles -------------------------------------------------------
enr = []
for c, d in df.groupby('condition'):
    d = d.assign(bin=pd.qcut(d.logm, 5, labels=['Q1 (lightest)', 'Q2', 'Q3', 'Q4', 'Q5 (heaviest)']))
    N = len(d)
    for k in CATS:
        K = int((d.category == k).sum())
        for b, db in d.groupby('bin', observed=True):
            n = len(db); x = int((db.category == k).sum())
            p_over = stats.hypergeom.sf(x - 1, N, K, n); p_under = stats.hypergeom.cdf(x, N, K, n)
            enr.append({'condition': c, 'category': k, 'mass_bin': b, 'observed': x, 'expected': K * n / N,
                        'obs/exp': x / (K * n / N) if K else np.nan, 'p_over': p_over, 'p_under': p_under})
enr = pd.DataFrame(enr)
enr['p_over_BH'] = bh(enr.p_over); enr['p_under_BH'] = bh(enr.p_under)
enr.round(4).to_csv(f'{OUT}/a02_M4c_mass_quintile_enrichment.csv', index=False)
pooled = enr.groupby(['category', 'mass_bin'], observed=True)[['observed', 'expected']].sum()
pooled['obs/exp'] = pooled.observed / pooled.expected
pooled['n_cond_over'] = enr.assign(s=enr.p_over_BH < 0.05).groupby(['category', 'mass_bin'], observed=True).s.sum()
pooled['n_cond_under'] = enr.assign(s=enr.p_under_BH < 0.05).groupby(['category', 'mass_bin'], observed=True).s.sum()
pooled.round(2).to_csv(f'{OUT}/a02_M4c_pooled.csv')
print('\nM4c - obs/exp by mass quintile, pooled over 11 conditions (n conditions significant over / under):')
print(pooled.round(2).to_string())

# ---- M4d: regression fitness ~ log10 mass x category + condition --------------------------------
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

reg = []
for k in CATS:
    d = df[df.category == k]
    if len(d) < 30:
        continue
    cond = pd.get_dummies(d.condition, drop_first=True, dtype=float)
    X = np.column_stack([np.ones(len(d)), d.logm.values, cond.values])
    beta, se, G = ols_cluster(d.fitness.values, X, d.bnum.values)
    t = beta[1] / se[1]
    within = d.assign(f=d.fitness - d.groupby('condition').fitness.transform('mean'),
                      m=d.logm - d.groupby('condition').logm.transform('mean'))
    rho = stats.spearmanr(within.m, within.f)
    reg.append({'category': k, 'n_obs': len(d), 'n_genes': G, 'slope_fitness_per_10x_mass': beta[1],
                'cluster_SE': se[1], 'p': 2 * stats.t.sf(abs(t), G - 1),
                'spearman_rho_within_condition': rho.statistic})
reg = pd.DataFrame(reg)
# all genes together: interaction model, mean-effect as reference
cat_d = pd.get_dummies(df.category, dtype=float)[['essential', 'important', 'enhancing']]
cond = pd.get_dummies(df.condition, drop_first=True, dtype=float)
X = np.column_stack([np.ones(len(df)), df.logm, cat_d.values, cat_d.values * df.logm.values[:, None], cond.values])
beta, se, G = ols_cluster(df.fitness.values, X, df.bnum.values)
inter = pd.DataFrame({'term': ['log10 mass (mean-effect slope)'] + [f'log10 mass x {k}' for k in cat_d.columns],
                      'estimate': [beta[1]] + list(beta[5:8]), 'cluster_SE': [se[1]] + list(se[5:8])})
inter['p'] = 2 * stats.t.sf(np.abs(inter.estimate / inter.cluster_SE), G - 1)
# continuous: all genes, no categories
Xc = np.column_stack([np.ones(len(df)), df.logm, cond.values])
bc, sc, _ = ols_cluster(df.fitness.values, Xc, df.bnum.values)
reg.to_csv(f'{OUT}/a02_M4d_regression_by_category.csv', index=False)
inter.to_csv(f'{OUT}/a02_M4d_interaction_model.csv', index=False)
print('\nM4d - fitness ~ log10 mass + condition, per category (cluster-robust SE by gene):')
print(reg.round(4).to_string())
print(inter.round(4).to_string())
print(f'All genes, no categories: slope = {bc[1]:.4f} per 10x mass (SE {sc[1]:.4f}), '
      f'Spearman within condition = {stats.spearmanr(df.logm - df.groupby("condition").logm.transform("mean"), df.fitness - df.groupby("condition").fitness.transform("mean")).statistic:.3f}')

# ---- m9: enhancing vs mean-effect, location and spread -----------------------------------------
m9 = []
for c, d in df.groupby('condition'):
    x = d.loc[d.category == 'enhancing', 'logm']; y = d.loc[d.category == 'mean-effect', 'logm']
    if len(x) < 5:
        continue
    m9.append({'condition': c, 'n_enhancing': len(x), 'median_enh': x.median(), 'median_mean_effect': y.median(),
               'IQR_enh': x.quantile(.75) - x.quantile(.25), 'IQR_mean_effect': y.quantile(.75) - y.quantile(.25),
               'SD_enh': x.std(), 'SD_mean_effect': y.std(),
               'p_location_MWU': stats.mannwhitneyu(x, y).pvalue,
               'p_spread_BrownForsythe': stats.levene(x, y, center='median').pvalue,
               'p_spread_AnsariBradley': stats.ansari(x - x.median(), y - y.median()).pvalue})
m9 = pd.DataFrame(m9)
for col in ['p_location_MWU', 'p_spread_BrownForsythe', 'p_spread_AnsariBradley']:
    m9[col + '_BH'] = bh(m9[col])
m9.round(4).to_csv(f'{OUT}/a02_m9_enhancing_vs_mean_effect.csv', index=False)
print('\nm9 - fitness-enhancing vs mean-effect (log10 mass):')
print(m9.round(3).to_string())

# ---- figure: mass distribution by category (pooled, log10) + obs/exp heat map -----------------
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={'width_ratios': [1, 1.3]})
for k in CATS:
    v = df.loc[df.category == k, 'logm']
    a1.hist(v, bins=60, density=True, histtype='step', lw=1.6, color=COLORS[k], label=f'{k} (n={len(v)})')
a1.set_xlabel('log10 protein mass (fg/cell)'); a1.set_ylabel('Density'); a1.legend(frameon=False, fontsize=8)
a1.set_title('Protein mass by fitness category (11 conditions pooled)')
hm = pooled['obs/exp'].unstack().reindex(CATS)
im = a2.imshow(np.log2(hm.values), cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto')
a2.set_xticks(range(hm.shape[1])); a2.set_xticklabels(hm.columns, rotation=30, ha='right')
a2.set_yticks(range(len(CATS))); a2.set_yticklabels(CATS)
for i in range(hm.shape[0]):
    for j in range(hm.shape[1]):
        a2.text(j, i, f'{hm.values[i, j]:.2f}', ha='center', va='center', fontsize=8)
plt.colorbar(im, ax=a2, label='log2 observed / expected')
a2.set_title('Category representation by protein-mass quintile')
fig.tight_layout(); fig.savefig(f'{OUT}/a02_mass_by_category.png', dpi=200)
