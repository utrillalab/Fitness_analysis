"""R1-M3 / R1-m3: ME-model predictions across all fitness categories; predicted flux vs measured mass.
R1-M8 : condition specificity of each fitness category (quantitative metric across conditions).
R1-m10: sensitivity of the fitness-enhancing gene set to the >=10-condition threshold, with
        enrichment of RegulonDB regulators; also without rpsG and rpsR, probable Tn-seq artefacts (R2-m2).
        Table S20p is built from a27_threshold_core.py (candidates and core)."""
import os, re
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from common import *

rng = np.random.default_rng(2)
df = long_table()
flux = load_flux()

# ---- M3 / m3 ------------------------------------------------------------------------------------
me = []
corr = []
for c in ME_CONDS:
    d = df[df.condition == c].copy()
    d['flux'] = d.bnum.map(flux[c])
    d['me_status'] = np.where(d.flux.isna(), 'absent from ME', np.where(d.flux == 0, 'zero flux', 'non-zero flux'))
    tot = d.mass_fg.sum()
    for k in CATS:
        dk = d[d.category == k]
        for s in ['absent from ME', 'zero flux', 'non-zero flux']:
            ds = dk[dk.me_status == s]
            me.append({'condition': c, 'category': k, 'me_status': s, 'n_genes': len(ds),
                       '%_of_category_genes': 100 * len(ds) / max(len(dk), 1),
                       '%_of_category_genes_in_ME': 100 * len(ds) / max((dk.me_status != 'absent from ME').sum(), 1)
                       if s != 'absent from ME' else np.nan,
                       'proteome_%': 100 * ds.mass_fg.sum() / tot})
    inme = d[d.me_status != 'absent from ME']
    # zero flux enriched among mean-effect genes? (Fisher, genes present in the ME-model)
    a = ((inme.category == 'mean-effect') & (inme.flux == 0)).sum(); b = ((inme.category == 'mean-effect') & (inme.flux > 0)).sum()
    cc = ((inme.category != 'mean-effect') & (inme.flux == 0)).sum(); dd = ((inme.category != 'mean-effect') & (inme.flux > 0)).sum()
    orr, pf = stats.fisher_exact([[a, b], [cc, dd]])
    nz = inme[(inme.flux > 0) & (inme.mass_fg > 0)]
    r_mass = stats.spearmanr(nz.flux, nz.mass_fg)
    r_fit = stats.spearmanr(inme.flux, inme.fitness)
    r_fit_nz = stats.spearmanr(nz.flux, nz.fitness)
    corr.append({'condition': c, 'genes_in_ME': len(inme), 'non_zero_flux_genes': len(nz),
                 'spearman_flux_vs_mass (non-zero)': r_mass.statistic, 'p_flux_mass': r_mass.pvalue,
                 'spearman_flux_vs_fitness (all ME genes)': r_fit.statistic, 'p_flux_fitness': r_fit.pvalue,
                 'spearman_flux_vs_fitness (non-zero)': r_fit_nz.statistic,
                 'OR_zero_flux_mean_effect_vs_other': orr, 'p_fisher': pf})
me = pd.DataFrame(me); corr = pd.DataFrame(corr)
me.round(2).to_csv(f'{OUT}/a03_M3_me_status_by_category.csv', index=False)
corr.round(4).to_csv(f'{OUT}/a03_m3_flux_correlations.csv', index=False)
print('M3 - % of genes present in the ME-model that have zero predicted flux, by category (median over 10 conditions):')
print(me[me.me_status == 'zero flux'].groupby('category')['%_of_category_genes_in_ME'].median().round(1).reindex(CATS).to_string())
print('\nM3 - proteome % by category and ME status (median over 10 conditions):')
print(me.groupby(['category', 'me_status']).proteome_pct.median().unstack().round(1).reindex(CATS).to_string()
      if 'proteome_pct' in me else me.groupby(['category', 'me_status'])['proteome_%'].median().unstack().round(1).reindex(CATS).to_string())
print('\nm3 / M3 correlations:')
print(corr.round(3).to_string())

# figure: ME status by category (median across conditions, % of category genes in ME)
piv = me[me.me_status != 'absent from ME'].groupby(['category', 'me_status'])['%_of_category_genes_in_ME'].median().unstack().reindex(CATS)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.5))
piv[['non-zero flux', 'zero flux']].plot.barh(stacked=True, ax=a1, color=['#b2182b', '#252525'])
a1.set_xlabel('% of genes present in the ME-model (median, 10 conditions)'); a1.set_ylabel('')
a1.set_title('Predicted use of proteins by fitness category'); a1.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2)
c = 'Glucose'; d = df[df.condition == c].assign(flux=lambda x: x.bnum.map(flux[c]))
d = d[(d.flux > 0) & (d.mass_fg > 0)]
for k in CATS:
    dk = d[d.category == k]; a2.scatter(np.log10(dk.flux), np.log10(dk.mass_fg), s=6, alpha=.6, color=COLORS[k], label=k)
r = stats.spearmanr(d.flux, d.mass_fg).statistic
a2.set_xlabel('log10 predicted translation flux (ME-model)'); a2.set_ylabel('log10 measured mass (fg/cell)')
a2.set_title(f'Glucose: predicted flux vs measured mass (Spearman ρ = {r:.2f})'); a2.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f'{OUT}/a03_me_model.png', dpi=200)

# ---- M8: condition specificity ------------------------------------------------------------------
fit, names = load_fitness()
cat11 = pd.DataFrame({c: classify(fit[c]) for c in CONDS}, index=fit.index)
spec = []
for k in CATS:
    sets = {c: set(cat11.index[cat11[c] == k]) for c in CONDS}
    jac = [len(sets[a] & sets[b]) / len(sets[a] | sets[b]) for i, a in enumerate(CONDS) for b in CONDS[i + 1:]
           if len(sets[a] | sets[b])]
    # null: same number of genes per condition, drawn at random
    nulls = []
    for _ in range(200):
        rs = {c: set(rng.choice(fit.index, len(sets[c]), replace=False)) for c in CONDS}
        nulls.append(np.mean([len(rs[a] & rs[b]) / len(rs[a] | rs[b]) for i, a in enumerate(CONDS) for b in CONDS[i + 1:]]))
    ever = (cat11 == k).sum(axis=1); ever = ever[ever > 0]
    spec.append({'category': k, 'mean_genes_per_condition': np.mean([len(s) for s in sets.values()]),
                 'genes_in_category_in_>=1_condition': len(ever),
                 'mean_pairwise_Jaccard_between_conditions': np.mean(jac), 'Jaccard_null_mean': np.mean(nulls),
                 'Jaccard_obs/null': np.mean(jac) / np.mean(nulls),
                 'median_conditions_per_gene (of 11)': ever.median(),
                 '%_genes_in_category_in_only_1_condition': 100 * (ever == 1).mean(),
                 '%_genes_in_category_in_all_11': 100 * (ever == 11).mean()})
spec = pd.DataFrame(spec).set_index('category')
spec.round(3).to_csv(f'{OUT}/a03_M8_condition_specificity.csv')
print('\nM8 - condition specificity per category (11 conditions):')
print(spec.round(3).to_string())
# pairwise test: conditions-per-gene distributions, important vs essential and vs mean-effect
cpg = {k: (cat11 == k).sum(axis=1)[lambda s: s > 0] for k in CATS}
for a, b in [('important', 'essential'), ('important', 'mean-effect'), ('important', 'enhancing')]:
    print(f'  conditions per gene, {a} vs {b}: MWU p = {stats.mannwhitneyu(cpg[a], cpg[b]).pvalue:.2e}, '
          f"Cliff's delta = {cliffs_delta(cpg[a], cpg[b]):.2f}")

# ---- m10: threshold for the recurrent fitness-enhancing set --------------------------------------
fa, names46 = load_fitness_all()
n_enh = (fa > 0.5).sum(axis=1)
net = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'NetworkRegulatorGene.tsv'),
                  sep='\t', comment='#', header=None)
reg_names = set(net[1].dropna().astype(str))
lower = {n.lower(): n for n in names46.dropna().unique()}
regulators = {lower[r.lower()] for r in reg_names if r.lower() in lower}
reg_b = set(names46.index[names46.isin(regulators)])
# Likely pooled-TnSeq artefacts (R2-m2): rpsG and rpsR, the only candidates essential in Keio, PEC and DEG (R1-m13)
artefact = set(names46.index[names46.isin(['rpsG', 'rpsR'])]) & set(fa.index)
print(f'\nRegulonDB regulators mapped to fitness genes: {len(reg_b)} of {len(names46)}; '
      f'probable artefacts (rpsG, rpsR) in fitness data: {len(artefact)}')
thr = []
N = len(fa); K = len(reg_b)
for k in range(3, 26):
    for label, universe in [('all genes', set(fa.index)), ('without rpsG and rpsR', set(fa.index) - artefact)]:
        sel = set(n_enh.index[n_enh >= k]) & universe
        x = len(sel & reg_b); Nu = len(universe); Ku = len(reg_b & universe)
        thr.append({'min_conditions': k, 'set': label, 'n_genes': len(sel), 'n_regulators': x,
                    '%_regulators': 100 * x / max(len(sel), 1), 'background_%': 100 * Ku / Nu,
                    'fold_enrichment': (x / max(len(sel), 1)) / (Ku / Nu),
                    'p_hypergeom': stats.hypergeom.sf(x - 1, Nu, Ku, len(sel)) if len(sel) else np.nan,
                    'flagged_artefact_genes': ', '.join(sorted(names46[list(set(n_enh.index[n_enh >= k]) & artefact)].astype(str)))
                    if label == 'all genes' else ''})
thr = pd.DataFrame(thr)
thr.round(4).to_csv(f'{OUT}/a03_m10_threshold_sensitivity.csv', index=False)
print(thr[thr.min_conditions.isin([5, 8, 10, 12, 15, 20])].drop(columns='flagged_artefact_genes').round(4).to_string())
print('Flagged genes in the >=10 set:', thr.loc[(thr.min_conditions == 10) & (thr.set == 'all genes'), 'flagged_artefact_genes'].iloc[0])

fig, ax = plt.subplots(figsize=(6.5, 4))
for label, ls in [('all genes', '-'), ('without rpsG and rpsR', '--')]:
    t = thr[thr.set == label]
    ax.plot(t.min_conditions, t.fold_enrichment, ls, color='#1b7837', label=f'fold enrichment ({label})')
ax2 = ax.twinx(); t = thr[thr.set == 'all genes']
ax2.bar(t.min_conditions, t.n_genes, color='#bdbdbd', alpha=.5); ax2.set_ylabel('genes in set (bars)')
ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
ax.axvline(10, color='k', lw=.8, ls=':'); ax.set_xlabel('minimum number of conditions (of 46) with fitness > 0.5')
ax.set_ylabel('regulator enrichment (fold)'); ax.legend(frameon=False, fontsize=7, loc='upper left')
fig.tight_layout(); fig.savefig(f'{OUT}/a03_m10_threshold.png', dpi=200)
