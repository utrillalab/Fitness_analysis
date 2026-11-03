"""R2-M2a / R1-M2: classification that accounts for per-gene measurement error (t-like statistic).

Source: Price et al. 2018 supplementary data for E. coli BW25113 (Keio), downloaded 2026-10-06 from
https://genomics.lbl.gov/supplemental/bigfit/html/Keio/ (fit_t.tab, fit_logratios_good.tab, expsUsed).
Replicate t scores are combined as in Price et al. 2018 Methods: t_comb = sum(t)/sqrt(n), divided by
sqrt((n^2 + n)/(2n)) when replicates share a start sample (conservative; used for all pairs here).
A gene keeps a non-mean-effect category only if |fitness| > 0.5 AND |t_comb| > 4 (Price et al. 2018 threshold).
Also: exact number of generations per experiment (expsUsed: Total.Generations) for the R1-M3 calculation.
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from common import *

P18 = os.path.join(os.path.dirname(HERE), 'Input_Data_Price2018')
tt = pd.read_csv(os.path.join(P18, 'fit_t.tab'), sep='\t').set_index('sysName')
lr = pd.read_csv(os.path.join(P18, 'fit_logratios_good.tab'), sep='\t').set_index('sysName')
ex = pd.read_csv(os.path.join(P18, 'expsUsed'), sep='\t')
ex['col'] = ex['name'] + ' ' + ex['short'].astype(str)

# replicate columns used for each condition (same as common.FIT_COL, plus the second replicate)
ft = pd.read_csv(os.path.join(REPO, 'Input_Data', 'fit_organism_Keio.tsv'), sep='\t').set_index('sysName')
REPS = {}
for c in CONDS:
    first = FIT_COL[c]
    if c == 'LB':
        REPS[c] = ['set2IT045 LB', 'set2IT069 LB']
    else:
        j = list(ft.columns).index(first)
        REPS[c] = [ft.columns[j], ft.columns[j + 1]]
print('replicates:', {c: [r.split()[0] for r in v] for c, v in REPS.items()})

# ---- 1. consistency with the repository data -----------------------------------------------------
chk = []
for c, reps in REPS.items():
    for r in reps:
        g = ft.index.intersection(lr.index)
        a = pd.to_numeric(ft.loc[g, r], errors='coerce'); b = pd.to_numeric(lr.loc[g, r], errors='coerce') if r in lr else None
        gen = ex.loc[ex.col == r, 'Total.Generations']
        chk.append({'condition': c, 'experiment': r, 'in_price2018_good': r in lr.columns,
                    'max_abs_diff_vs_repo': float(np.nanmax(np.abs(a - b))) if b is not None else np.nan,
                    'pearson_vs_repo': float(np.corrcoef(a[b.notna()], b[b.notna()])[0, 1]) if b is not None else np.nan,
                    'generations': float(gen.iloc[0]) if len(gen) else np.nan,
                    'start_OD': float(ex.loc[ex.col == r, 'StartOD'].iloc[0]) if len(gen) else np.nan,
                    'end_OD': float(ex.loc[ex.col == r, 'EndOD'].iloc[0]) if len(gen) else np.nan})
chk = pd.DataFrame(chk)
chk.round(4).to_csv(f'{OUT}/a05_consistency_and_generations.csv', index=False)
print(chk.round(3).to_string())

# ---- 2. combined t and t-aware classification ---------------------------------------------------
fit, names = load_fitness()
tcomb = pd.DataFrame({c: tt[REPS[c]].sum(axis=1) / np.sqrt(2) / np.sqrt(1.5) for c in CONDS})
tcomb = tcomb.reindex(fit.index)


def classify_t(f, t, tmin=4.0):
    base = classify(f)
    sig = np.abs(np.asarray(t, float)) > tmin
    return np.where((base != 'mean-effect') & (base != 'NA') & ~sig, 'mean-effect', base)


df = long_table()
df['t'] = [tcomb.at[b, c] for b, c in zip(df.bnum, df.condition)]
df['cat_t'] = classify_t(df.fitness, df.t)
df['sig'] = np.abs(df.t) > 4

# reclassification counts (all 3,789 genes, not only those with proteomics)
rc = []
for c in CONDS:
    base = pd.Series(classify(fit[c]), index=fit.index); new = pd.Series(classify_t(fit[c], tcomb[c]), index=fit.index)
    small_sig = ((base == 'mean-effect') & (np.abs(tcomb[c]) > 4)).sum()
    rc.append({'condition': c, **{f'{k}_before': int((base == k).sum()) for k in CATS},
               **{f'{k}_after': int((new == k).sum()) for k in CATS},
               'moved_to_mean_effect': int(((base != 'mean-effect') & (new == 'mean-effect')).sum()),
               'mean_effect_with_|t|>4 (small but significant)': int(small_sig)})
rc = pd.DataFrame(rc).set_index('condition')
rc.to_csv(f'{OUT}/a05_reclassification_counts.csv')
print('\nGenes per category before -> after requiring |t| > 4:')
print(rc.to_string())

# proteome fractions: original vs t-aware
def fr(d, col):
    tot = d.groupby('condition').mass_fg.sum()
    return (d.groupby(['condition', col]).mass_fg.sum().unstack().div(tot, axis=0) * 100).reindex(columns=CATS).fillna(0)

frac = pd.concat({'cut-offs only': fr(df, 'category'), 'cut-offs + |t| > 4': fr(df, 'cat_t')}, axis=1).round(1)
frac.to_csv(f'{OUT}/a05_proteome_fractions_t.csv')
print('\nProteome fraction (%) by category:')
print(frac.to_string())
me = df[df.category == 'mean-effect']
ms = me.groupby('condition').apply(lambda d: 100 * d.loc[d.sig, 'mass_fg'].sum() / df[df.condition == d.name].mass_fg.sum())
print('\nMean-effect genes with |t| > 4 (small but significant effect), % of proteome:', ms.round(1).to_dict())

# ---- 3. key results re-tested with the t-aware classification ----------------------------------
# 3a mass by category (Cliff's delta vs mean-effect)
d2 = df[df.mass_fg > 0].assign(logm=lambda x: np.log10(x.mass_fg))
eff = []
for c, d in d2.groupby('condition'):
    ref = d.loc[d.cat_t == 'mean-effect', 'logm']
    for k in ['essential', 'important', 'enhancing']:
        x = d.loc[d.cat_t == k, 'logm']
        if len(x) >= 3:
            eff.append({'condition': c, 'category': k, 'n': len(x), 'cliffs_delta': cliffs_delta(x, ref),
                        'p': stats.mannwhitneyu(x, ref).pvalue})
eff = pd.DataFrame(eff); eff['p_holm'] = holm(eff.p)
eff.round(4).to_csv(f'{OUT}/a05_effect_sizes_t.csv', index=False)
print('\nCliff\'s delta vs mean-effect with t-aware categories (median; n conditions Holm<0.05):')
print(eff.groupby('category').agg(delta=('cliffs_delta', 'median'), sig=('p_holm', lambda p: int((p < 0.05).sum())),
                                  n=('p_holm', 'size')).round(2).to_string())

# 3b recurrent fitness-enhancing genes over 46 conditions with |t| > 4, and regulator enrichment
fa, names46 = load_fitness_all()
cols46 = list(fa.columns)
t46 = {}
for col in cols46:
    if col == 'LB':
        t46[col] = tt[REPS['LB']].sum(axis=1) / np.sqrt(3)
    else:
        j = list(ft.columns).index(col)
        pair = [ft.columns[j], ft.columns[j + 1]]
        t46[col] = tt[pair].sum(axis=1) / np.sqrt(3) if all(p in tt.columns for p in pair) else np.nan
t46 = pd.DataFrame(t46).reindex(fa.index)
n_enh = (fa > 0.5).sum(axis=1)
n_enh_t = ((fa > 0.5) & (t46 > 4)).sum(axis=1)
net = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'NetworkRegulatorGene.tsv'), sep='\t', comment='#', header=None)
lower = {n.lower(): n for n in names46.dropna().unique()}
regs = {lower[r.lower()] for r in set(net[1].dropna().astype(str)) if r.lower() in lower}
reg_b = set(names46.index[names46.isin(regs)])
N, K = len(fa), len(reg_b)
out = []
for label, n in [('fitness > 0.5', n_enh), ('fitness > 0.5 and t > 4', n_enh_t)]:
    for k in [5, 10, 15]:
        sel = set(n.index[n >= k]); x = len(sel & reg_b)
        out.append({'criterion': label, 'min_conditions': k, 'n_genes': len(sel), 'n_regulators': x,
                    'fold_enrichment': (x / max(len(sel), 1)) / (K / N),
                    'p_hypergeom': stats.hypergeom.sf(x - 1, N, K, len(sel)) if sel else np.nan})
out = pd.DataFrame(out)
out.round(4).to_csv(f'{OUT}/a05_enhancing_regulators_t.csv', index=False)
print('\nRecurrent fitness-enhancing genes and regulator enrichment:')
print(out.round(4).to_string())
kept = set(n_enh_t.index[n_enh_t >= 10]); orig = set(n_enh.index[n_enh >= 10])
print(f'Of the 108 genes (>=10 conditions), {len(orig & kept)} remain with t > 4 in >=10 conditions; '
      f'rpoS: {n_enh.get("b2741")} -> {n_enh_t.get("b2741")} conditions')

# ---- 4. figure: fitness vs t, one representative condition + pooled counts ----------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
c = 'Glucose'
f = fit[c]; t = tcomb[c]
ax = axes[0]
sig = np.abs(t) > 4
ax.scatter(f[~sig], t[~sig], s=4, color='#bdbdbd', label='|t| ≤ 4')
ax.scatter(f[sig], t[sig], s=4, color='#2166ac', label='|t| > 4')
for x in (-0.5, 0.5): ax.axvline(x, color='k', ls='--', lw=0.8)
ax.axvline(-2, color='k', ls=':', lw=0.8)
for y in (-4, 4): ax.axhline(y, color='#2166ac', ls='-', lw=0.6)
ax.set_xlim(-4, 2); ax.set_ylim(-30, 15)
ax.set_xlabel('Fitness (glucose)'); ax.set_ylabel('Combined t-like statistic'); ax.legend(frameon=False, fontsize=8)
ax.set_title('Glucose: fitness vs t (Price et al. 2018 data)')
ax = axes[1]
fr2 = frac.xs('mean-effect', axis=1, level=1)
x = np.arange(len(CONDS)); w = 0.38
ax.bar(x - w / 2, fr2.loc[CONDS, 'cut-offs only'], w, color='#636363', label='cut-offs only')
ax.bar(x + w / 2, fr2.loc[CONDS, 'cut-offs + |t| > 4'], w, color='#2166ac', label='cut-offs + |t| > 4')
ax.set_xticks(x); ax.set_xticklabels(CONDS, rotation=50, ha='right')
ax.set_ylabel('Mean-effect proteome fraction (%)'); ax.legend(frameon=False, fontsize=8)
ax.set_title('Mean-effect fraction with and without the t filter')
for a in axes: a.spines[['top', 'right']].set_visible(False)
fig.tight_layout(); fig.savefig(f'{OUT}/a05_t_statistic.png', dpi=200)

# ---- 5. candidate (108) vs significant core recurrent fitness-enhancing genes
# Core (option B, chosen 2026-10-07): recurrent (fitness > 0.5 in >= 10 conditions) AND significant (t > 4) in >= 1 of them ----------------
eco = pd.read_csv(os.path.join(REPO, 'Input_Data', 'All-genes-of-E.-coli-K-12-substr-MG1655.txt'), sep='\t')
prod = eco.drop_duplicates('Accession-1').set_index('Accession-1').Product
transl = set(eco.loc[eco.Product.astype(str).str.contains('ribosomal subunit protein') |
                     eco.Product.astype(str).str.match(r'^tRNA-[A-Z][a-z]{2}\('), 'Accession-1'])
cand = pd.DataFrame({'bnumber': sorted(orig)})
cand['gene'] = cand.bnumber.map(names46)
cand['product'] = cand.bnumber.map(prod)
cand['n_conditions_fitness>0.5 (of 46)'] = cand.bnumber.map(n_enh)
cand['n_conditions_fitness>0.5_and_t>4'] = cand.bnumber.map(n_enh_t)
cand['significant_core'] = cand['n_conditions_fitness>0.5_and_t>4'] >= 1
cand['RegulonDB_regulator'] = cand.bnumber.isin(reg_b)
cand['ribosomal_protein_or_tRNA'] = cand.bnumber.isin(transl)
cand = cand.sort_values(['significant_core', 'n_conditions_fitness>0.5_and_t>4', 'n_conditions_fitness>0.5 (of 46)'],
                        ascending=False)
cand.to_csv(f'{OUT}/a05_enhancing_108_candidates_32_core.csv', index=False)
core = cand[cand.significant_core]
print(f'\nSignificant core: {len(core)} genes; regulators {int(core.RegulonDB_regulator.sum())}; '
      f'r-protein/tRNA {int(core.ribosomal_protein_or_tRNA.sum())} (of {int(cand.ribosomal_protein_or_tRNA.sum())} among the 108)')
print(core[['gene', 'n_conditions_fitness>0.5 (of 46)', 'n_conditions_fitness>0.5_and_t>4', 'RegulonDB_regulator', 'product']].to_string(index=False))
print('\nr-protein/tRNA candidates:', cand.loc[cand.ribosomal_protein_or_tRNA, ['gene', 'n_conditions_fitness>0.5 (of 46)', 'n_conditions_fitness>0.5_and_t>4']].to_string(index=False))
