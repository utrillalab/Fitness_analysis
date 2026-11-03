"""Fig. 6 redrawn at 300 dpi (R2-M2b): fitness-enhancing regulators.

(a, b) Proteome percentage of the RpoS (sigma38) regulon, all 340 targets and the 131 targets unique to RpoS, per carbon
source: same calculation as Functional_Analysis/Fig06_rpoS_sum_pf.ipynb (RegulonDB sigma network; proteome fractions
relative to the summed mass of all genes in the sigma network, as stated in Methods).
(c) Number of TFs regulating each gene in the whole TRN and in the genes regulated by the five fitness-enhancing TFs
(RbsR, GcvA, Fis, SlyA, RcsB), from D. Ledezma-Tejeida's RegulonDB 11.2 network (R1_analyses/input_DLT/NetWorkTFGene.txt);
replaces the former panel c, which also showed a single random sample of 323 genes. Two-sided Wilcoxon-Mann-Whitney test.
Outputs (R1_analyses/output/): Fig06_revised.{png,svg,pdf}; fig6_values.csv
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('FA_REPO', os.path.dirname(HERE))
OUT = os.path.join(HERE, 'R1_analyses', 'output')
INF = os.path.join(REPO, 'Functional_Analysis', 'Infiles')
from npj_style import *   # npj style (28 Oct 2026): Arial 7 pt at print size, white background
CONDS = ['Glucose', 'LB', 'Glucosamine', 'Glycerol', 'Xylose', 'Mannose', 'Galactose', 'Fructose', 'Acetate',
         'Pyruvate', 'Succinate']

# --- a, b: RpoS regulon (as in the notebook) ---
sig = pd.read_csv(os.path.join(INF, 'NetworkSigmaGene.tsv'), sep='\t', comment='#')
sig = sig.rename(columns={'1)sigmaName': 'Sigma', '2)regulatedName': 'Target'})[['Sigma', 'Target']]
prot = pd.read_csv(os.path.join(INF, 'allgen4691_fitfg_newname.csv'))[['gene'] + CONDS].rename(columns={'gene': 'Target'})
prot = prot[prot.Target.notnull()].drop_duplicates('Target')
net = pd.merge(sig, prot, on='Target', how='left').fillna(0)
tot = net[CONDS].sum()
s38 = net[net.Sigma == 'sigma38']
others = set(net.loc[net.Sigma != 'sigma38', 'Target'])
uniq = net[(net.Sigma == 'sigma38') & ~net.Target.isin(others)]
assert s38.Target.nunique() == 340 and uniq.Target.nunique() == 131, (s38.Target.nunique(), uniq.Target.nunique())
pa = (100 * s38[CONDS].sum() / tot).sort_values(ascending=False)
pb = (100 * uniq[CONDS].sum() / tot).sort_values(ascending=False)

# --- c: TFs per gene (D. Ledezma-Tejeida, RegulonDB 11.2) ---
n = pd.read_csv(os.path.join(HERE, 'R1_analyses', 'input_DLT', 'NetWorkTFGene.txt'), sep='\t', comment='#', header=0, dtype=str)
n.columns = range(n.shape[1])
E = n[[1, 4]].apply(lambda s: s.str.strip()).drop_duplicates().rename(columns={1: 'tf', 4: 'gene'})
deg = E.groupby('gene').tf.nunique()
FE = ['RbsR', 'GcvA', 'Fis', 'SlyA', 'RcsB']
T = sorted(set(E[E.tf.isin(FE)].gene))
d_all, d_fe = deg.values, deg.reindex(T).values
p = stats.mannwhitneyu(d_fe, d_all, alternative='two-sided').pvalue
print(f'TRN {len(d_all)} genes, mean {d_all.mean():.2f}, median {np.median(d_all):g}, >=3 {100*(d_all>=3).mean():.1f}%')
print(f'FE  {len(d_fe)} genes, mean {d_fe.mean():.2f}, median {np.median(d_fe):g}, >=3 {100*(d_fe>=3).mean():.1f}%;  WMW p = {p:.1e}')
for name, v in pa.items():
    print(f'a {name:12s} {v:5.2f}   b {pb.index[list(pa.index).index(name)]:12s} {pb.iloc[list(pa.index).index(name)]:5.2f}')

fig = plt.figure(figsize=(DOUBLE, 120 * MM))
gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], hspace=0.95, wspace=0.3, left=0.08, right=0.98, bottom=0.17, top=0.94)
BLUE = PAL['enhancing']   # RpoS is a fitness-enhancing regulator
axs = []
for k, (ser, title) in enumerate([(pa, 'All RpoS-dependent genes (n = 340)'), (pb, 'Genes unique to RpoS (n = 131)')]):
    ax = fig.add_subplot(gs[k, 0]); axs.append(ax)
    ax.bar(range(len(ser)), ser.values, color=BLUE, width=0.75)
    ax.set_xticks(range(len(ser))); ax.set_xticklabels(ser.index, rotation=55, ha='right', rotation_mode='anchor')
    ax.set_ylabel('Proteome (%)'); ax.set_title(title)
ax = fig.add_subplot(gs[:, 1]); axs.append(ax)
rng = np.random.default_rng(6)
# points in grey (all TRN genes) and blue (targets of the fitness-enhancing TFs): no red-green pair
for i, (vals, col) in enumerate([(d_all, PAL['mean-effect']), (d_fe, PAL['enhancing'])]):
    ax.scatter(i + rng.uniform(-0.28, 0.28, len(vals)), vals + rng.uniform(-0.25, 0.25, len(vals)), s=1.5, alpha=0.4,
               color=col, lw=0, rasterized=True)
ax.boxplot([d_all, d_fe], positions=[0, 1], widths=0.55, showfliers=False, patch_artist=False, whis=1.5,
           medianprops=dict(color='black', lw=1.2), boxprops=dict(color='#333333', lw=0.8),
           whiskerprops=dict(color='#333333', lw=0.8), capprops=dict(color='#333333', lw=0.8))
ax.set_xticks([0, 1]); ax.set_xticklabels([f'All TRN genes\n(n = {len(d_all):,})', f'Genes regulated by\nfitness-enhancing TFs\n(n = {len(d_fe)})'])
ax.set_ylabel('Number of TFs regulating each gene'); ax.set_ylim(0, max(d_all.max(), d_fe.max()) + 2)
y = max(d_all.max(), d_fe.max()) + 0.8
ax.plot([0, 0, 1, 1], [y - 0.3, y, y, y - 0.3], color='black', lw=0.6)
ax.text(0.5, y + 0.1, '***', ha='center', va='bottom')
for ax, lt in zip(axs, 'abc'):
    letter(fig, ax, lt, dx=-0.01)
save(fig, os.path.join(OUT, 'Fig06_revised'))
pd.DataFrame({'condition': pa.index, 'all_340_pct': pa.values}).merge(
    pd.DataFrame({'condition': pb.index, 'unique_131_pct': pb.values}), on='condition').to_csv(
    os.path.join(OUT, 'fig6_values.csv'), index=False)
print('saved Fig06_revised')
