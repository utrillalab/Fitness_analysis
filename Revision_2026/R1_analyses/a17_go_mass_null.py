"""R1-M5b: does the proteome mass of each enriched GO term exceed what its number of genes predicts?

For every term enriched in a fitness category (corrected ViSEAGO results, output/viseago/corrected/), the
observed mass fraction = summed mass of the term's genes in that category / total measured proteome.
Null: the same number of genes (those with proteomics data) drawn at random 10,000 times from
  (a) all measured proteins in that condition, and
  (b) the measured proteins of the same fitness category.
Null sums for every set size come from one cumulative sum per permutation. Masses are not divided by the
number of GO annotations (R1-m7)."""
import os, re
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *

NPERM = 10000
rng = np.random.default_rng(17)
VIS = os.path.join(OUT, 'viseago', 'corrected')
CAT_COL = {'Essential': 'essential', 'Important': 'important', 'No_effect': 'mean-effect',
           'Fitness_enhancing': 'enhancing'}

fit, names = load_fitness()
mass = load_mass()
df = long_table()
eco = pd.read_csv(os.path.join(REPO, 'Input_Data', 'All-genes-of-E.-coli-K-12-substr-MG1655.txt'), sep='\t')
sym2b = {}
for n, b in zip(eco['Gene Name'], eco['Accession-1']):
    if isinstance(n, str) and isinstance(b, str):
        sym2b.setdefault(n.lower(), b)
for b, n in names.items():
    if isinstance(n, str):
        sym2b.setdefault(n.lower(), b)


def null_sums(m):
    """(NPERM x len(m)) cumulative sums of randomly permuted masses: column k-1 = sum of k random genes."""
    out = np.empty((NPERM, len(m)))
    for i in range(NPERM):
        out[i] = np.cumsum(rng.permutation(m))
    return out


rows, unmatched = [], set()
for c in CONDS:
    f = os.path.join(VIS, f'BP_Results_{c}.xls')
    if not os.path.exists(f):
        print('missing', f); continue
    res = pd.read_csv(f, sep='\t')
    m_all = mass[c].dropna()
    tot = m_all.sum()
    d = df[df.condition == c].set_index('bnum')
    null_all = null_sums(m_all.values)
    null_cat = {}
    for col, cat in CAT_COL.items():
        genes_cat = d.index[(d.category == cat) & d.mass_fg.notna()]
        null_cat[cat] = (null_sums(d.loc[genes_cat, 'mass_fg'].values), d.loc[genes_cat, 'mass_fg'].sum())
        sig = res[res[f'{col}.Significant_genes_symbol'].notna()]
        for _, r in sig.iterrows():
            syms = [s for s in str(r[f'{col}.Significant_genes_symbol']).split(';') if s]
            bs = [sym2b.get(s.lower()) for s in syms]
            unmatched |= {s for s, b in zip(syms, bs) if b is None}
            with_mass = sorted(b for b in set(bs) if b in m_all.index)  # sorted: float sums depend on the order
            k = len(with_mass)
            if k == 0:
                continue
            obs = m_all.loc[with_mass].sum()
            na = null_all[:, k - 1]
            ncat, cat_tot = null_cat[cat]
            row = {'condition': c, 'category': cat, 'GO.ID': r['GO.ID'], 'term': r['term'],
                   'n_genes_term_category': len(syms), 'n_with_proteomics': k,
                   'observed_mass_%': 100 * obs / tot, 'null_all_mean_%': 100 * na.mean() / tot,
                   'obs/exp_vs_all_proteins': obs / na.mean(),
                   'p_vs_all (upper)': (1 + (na >= obs).sum()) / (NPERM + 1),
                   'p_vs_all (lower)': (1 + (na <= obs).sum()) / (NPERM + 1)}
            if k <= ncat.shape[1]:
                nc = ncat[:, k - 1]
                row.update({'obs/exp_vs_same_category': obs / nc.mean(),
                            'p_vs_category (upper)': (1 + (nc >= obs).sum()) / (NPERM + 1),
                            'p_vs_category (lower)': (1 + (nc <= obs).sum()) / (NPERM + 1)})
            rows.append(row)
    print(c, 'done')
res = pd.DataFrame(rows)
for side in ['upper', 'lower']:
    res[f'p_vs_all ({side})_BH'] = bh(res[f'p_vs_all ({side})'])
    ok = res[f'p_vs_category ({side})'].notna()
    res.loc[ok, f'p_vs_category ({side})_BH'] = bh(res.loc[ok, f'p_vs_category ({side})'])
res.to_csv(f'{OUT}/a17_go_mass_null.csv', index=False)
if unmatched:
    print('symbols not matched (excluded):', len(unmatched), sorted(unmatched)[:20])

# summary per category
def summ(g):
    return pd.Series({'terms': len(g),
                      'median obs/exp vs all': g['obs/exp_vs_all_proteins'].median(),
                      'heavier than random (BH<0.05)': int((g['p_vs_all (upper)_BH'] < 0.05).sum()),
                      'lighter than random (BH<0.05)': int((g['p_vs_all (lower)_BH'] < 0.05).sum()),
                      'median obs/exp vs same category': g['obs/exp_vs_same_category'].median(),
                      'heavier than category (BH<0.05)': int((g['p_vs_category (upper)_BH'] < 0.05).sum()),
                      'lighter than category (BH<0.05)': int((g['p_vs_category (lower)_BH'] < 0.05).sum())})
s = res.groupby('category').apply(summ).reindex(CATS)
s.to_csv(f'{OUT}/a17_go_mass_null_summary.csv')
print('\nAll conditions, all enriched terms:'); print(s.round(2).to_string())

# Fig. 4 terms: top 5 by mass per category (top 2 for enhancing), glucose and galactose
top = []
for c in ['Glucose', 'Galactose']:
    for cat in CATS:
        g = res[(res.condition == c) & (res.category == cat)].sort_values('observed_mass_%', ascending=False)
        top.append(g.head(2 if cat == 'enhancing' else 5))
top = pd.concat(top)
top.to_csv(f'{OUT}/a17_fig4_terms_null.csv', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 45)
print('\nFig. 4 terms (top by mass):')
print(top[['condition', 'category', 'term', 'n_with_proteomics', 'observed_mass_%', 'obs/exp_vs_all_proteins',
           'p_vs_all (upper)_BH', 'p_vs_all (lower)_BH', 'obs/exp_vs_same_category']].round(3).to_string(index=False))

# figure: observed vs expected mass for every enriched term (log-log), coloured by category
PAL = {'essential': '#b2182b', 'important': '#e08214', 'mean-effect': '#969696', 'enhancing': '#2166ac'}
fig, ax = plt.subplots(figsize=(6, 5.5))
for cat in CATS:
    g = res[res.category == cat]
    ax.scatter(g['null_all_mean_%'], g['observed_mass_%'], s=10, alpha=0.6, color=PAL[cat], label=f'{cat} ({len(g)})')
lim = [res[['null_all_mean_%', 'observed_mass_%']].min().min() * 0.8, res[['null_all_mean_%', 'observed_mass_%']].max().max() * 1.2]
ax.plot(lim, lim, 'k--', lw=0.8)
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel('Expected mass, random genes of the same number (% of proteome)')
ax.set_ylabel('Observed mass of the term (% of proteome)')
ax.set_title('Enriched GO terms, 11 conditions'); ax.legend(frameon=False, fontsize=8)
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()
for ext in ('png', 'svg'):
    fig.savefig(f'{OUT}/a17_go_mass_null.{ext}', dpi=300 if ext == 'png' else None)
