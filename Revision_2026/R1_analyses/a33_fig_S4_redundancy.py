"""Supplementary Figure S4 (TF redundancy, R2-M2b), redrawn from D. Ledezma-Tejeida's network and null models.

Network: RegulonDB 11.2 TF-gene network (input_DLT/NetWorkTFGene.txt); without global TFs: input_DLT/TRN_noglobals.txt.
Fitness-enhancing (FE) TFs: RbsR, GcvA, Fis, SlyA, RcsB. Statistic: mean number of TFs regulating each target.
(A) Distribution of the number of TFs per gene in the whole TRN and in the genes regulated by the FE TFs (% of genes).
(B) Null distributions of the statistic: 10,000 uniform random sets of the same number of genes, and 10,000
    degree-preserving randomizations of the network (curveball; as in her random_networks script: seed 31, burn-in 20,000
    trades, 1,000 trades between saved networks); dashed line, observed value.
(C) TFs per gene in the TRN and in the FE targets after removing the seven global TFs (ArcA, CRP, Fis, FNR, H-NS, IHF,
    Lrp); two-sided Wilcoxon-Mann-Whitney test.
Outputs (output/): Fig_S4_redundancy.{png,svg,pdf}; a33_nulls.csv (null values); a33_summary.csv
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'output')
IN = os.path.join(HERE, 'input_DLT')
plt.rcParams.update({'svg.fonttype': 'none', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
FE = ['RbsR', 'GcvA', 'Fis', 'SlyA', 'RcsB']
B, TRADES, BURNIN = 10000, 1000, 20000
rng = np.random.default_rng(31)


def load(f):
    n = pd.read_csv(os.path.join(IN, f), sep='\t', comment='#', header=0, dtype=str)
    n.columns = range(n.shape[1])
    return n[[1, 4]].apply(lambda s: s.str.strip()).drop_duplicates().rename(columns={1: 'tf', 4: 'gene'})


def curveball(adj, trades):
    k = len(adj)
    for _ in range(trades):
        i, j = rng.choice(k, 2, replace=False)
        a, b = adj[i], adj[j]
        ai, bj = sorted(a - b), sorted(b - a)  # sorted: set order depends on the hash seed of each run
        if not ai or not bj:
            continue
        pool = ai + bj
        rng.shuffle(pool)
        common = a & b
        adj[i] = common | set(pool[:len(ai)])
        adj[j] = common | set(pool[len(ai):])


E = load('NetWorkTFGene.txt')
deg = E.groupby('gene').tf.nunique()
targ = E.groupby('tf').gene.apply(set)
T = set().union(*[targ[t] for t in FE])
obs = deg.reindex(list(T)).mean()
NULLS = os.path.join(OUT, 'a33_nulls.csv')
if os.path.exists(NULLS) and '--recompute' not in os.sys.argv:  # reuse the saved nulls (10,000 networks take ~10 min)
    _n = pd.read_csv(NULLS); uni, net = _n.uniform.values, _n.degree_preserving.values
else:
    uni = np.array([rng.choice(deg.values, len(T), replace=False).mean() for _ in range(B)])
    tfs = list(targ.index)
    adj = [set(targ[t]) for t in tfs]
    idx = [tfs.index(t) for t in FE]
    curveball(adj, BURNIN)
    net = []
    for b in range(B):
        curveball(adj, TRADES)
        net.append(deg.reindex(list(set().union(*[adj[i] for i in idx]))).mean())
        if b % 1000 == 0:
            print(b, flush=True)
    net = np.array(net)
p_uni = (1 + (uni >= obs).sum()) / (B + 1)
p_net = (1 + (net >= obs).sum()) / (B + 1)

Eg = load('TRN_noglobals.txt')
dg = Eg.groupby('gene').tf.nunique()
Tg = set(Eg[Eg.tf.isin(FE)].gene)
pg = stats.mannwhitneyu(dg.reindex(sorted(Tg)).values, dg.values, alternative='two-sided').pvalue
summ = pd.DataFrame([{'observed_mean': obs, 'n_targets': len(T), 'n_genes': len(deg),
                      'uniform_mean': uni.mean(), 'uniform_max': uni.max(), 'p_uniform': p_uni,
                      'network_mean': net.mean(), 'network_max': net.max(), 'p_network': p_net,
                      'noglobal_TRN_genes': len(dg), 'noglobal_TRN_mean': dg.mean(), 'noglobal_FE_genes': len(Tg),
                      'noglobal_FE_mean': dg.reindex(sorted(Tg)).mean(), 'noglobal_WMW_p': pg}])
summ.to_csv(os.path.join(OUT, 'a33_summary.csv'), index=False)
pd.DataFrame({'uniform': uni, 'degree_preserving': net}).to_csv(os.path.join(OUT, 'a33_nulls.csv'), index=False)
print(summ.T.to_string())

G, O = '#66c2a5', '#fc8d62'
fig, axs = plt.subplots(1, 3, figsize=(10, 3.4), gridspec_kw={'width_ratios': [1.25, 1.25, 1], 'wspace': 0.35})
ax = axs[0]
bins = np.arange(0.5, 17.5, 1)
for vals, col, lab, off in [(deg.values, G, f'TRN (n = {len(deg):,})', -0.2),
                            (deg.reindex(list(T)).values, O, 'targets of fitness-\nenhancing TFs ' + f'(n = {len(T)})', 0.2)]:
    h, _ = np.histogram(vals, bins)
    ax.bar(np.arange(1, 17) + off, 100 * h / h.sum(), width=0.4, color=col, label=lab)
ax.set_xlabel('Number of TFs regulating each gene'); ax.set_ylabel('Genes (%)'); ax.set_xticks(range(1, 17, 2))
ax.set_ylim(0, 50); ax.legend(frameon=False, fontsize=7.5, loc='lower right', bbox_to_anchor=(1.0, 0.62))
ax = axs[1]
ax.hist(uni, bins=50, color='#999999', alpha=0.8, label='random gene sets')
ax.hist(net, bins=50, color='#8da0cb', alpha=0.9, label='degree-preserving networks')
ax.axvline(obs, color='black', ls='--', lw=1)
ax.text(obs, ax.get_ylim()[1] * 0.95, f' observed\n {obs:.2f}', fontsize=7.5, va='top')
ax.set_xlabel('Mean number of TFs per target'); ax.set_ylabel(f'Null replicates (of {B:,})')
ax.legend(frameon=False, fontsize=7.5, loc='upper left')
ax = axs[2]
r = np.random.default_rng(6)
for i, (vals, col) in enumerate([(dg.values, G), (dg.reindex(sorted(Tg)).values, O)]):
    ax.scatter(i + r.uniform(-0.28, 0.28, len(vals)), vals + r.uniform(-0.25, 0.25, len(vals)), s=3, alpha=0.4,
               color=col, lw=0, rasterized=True)
ax.boxplot([dg.values, dg.reindex(sorted(Tg)).values], positions=[0, 1], widths=0.55, showfliers=False,
           medianprops=dict(color='#333333', lw=1.5), boxprops=dict(color='#555555'),
           whiskerprops=dict(color='#555555'), capprops=dict(color='#555555'))
ax.set_xticks([0, 1]); ax.set_xticklabels([f'TRN\n(n = {len(dg):,})', f'fitness-enhancing\nTF targets (n = {len(Tg)})'],
                                          fontsize=8)
ax.set_ylabel('TFs per gene, without global TFs')
for a, l in zip(axs, 'abc'):  # npj: lower-case bold
    a.text(-0.18, 1.04, l, transform=a.transAxes, fontsize=12, fontweight='bold')
for ext in ('png', 'svg', 'pdf'):
    fig.savefig(os.path.join(OUT, f'Fig_S4_redundancy.{ext}'), dpi=300, bbox_inches='tight')
print('saved')
