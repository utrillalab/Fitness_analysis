"""R2-M2b: null models for the regulatory redundancy of the targets of the five fitness-enhancing TFs.

Network: TF-gene interactions of RegulonDB in the repository (Functional_Analysis/Infiles/NetworkRegulatorGene.tsv, rows
with a TF regulator, RDBECOLITFC*); in-degree of a gene = number of TFs regulating it. Fitness-enhancing (FE) TFs:
RbsR, GcvA, Fis, SlyA, RcsB. Statistics of their target set: mean in-degree and % of targets with >= 3 TFs.
Nulls (empirical one-sided p = (1 + #null >= observed) / (n + 1)):
  1. uniform: random gene sets of the same size (the Monte Carlo of D. Ledezma-Tejeida's text);
  2. matched TFs: 5 random other TFs, each matched to one FE TF by number of targets (the 6 nearest in out-degree);
  3. degree-preserving network: curveball randomization of the bipartite TF-gene network, keeping the out-degree of every
     TF and the in-degree of every gene (the igraph-type null the reviewer asked for); statistic recomputed on the
     targets of the five FE TFs in each randomized network.
Leave-one-out: (a) each FE TF removed from the FE set; (b) each global TF (ArcA, CRP, Fis, FNR, H-NS, IHF, Lrp) removed
from the whole network, and all seven together; null 3 recomputed each time.
Output: output/a31_tf_redundancy_null.csv
"""
import os
import numpy as np, pandas as pd
from common import REPO, OUT

rng = np.random.default_rng(31)
N_UNI, N_TF, N_NET, TRADES = 10000, 10000, 1000, 1000
FE = ['RbsR', 'GcvA', 'Fis', 'SlyA', 'RcsB']
GLOBAL = ['ArcA', 'CRP', 'Fis', 'FNR', 'H-NS', 'IHF', 'Lrp']

n = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'NetworkRegulatorGene.tsv'), sep='\t',
                comment='#', header=None)
n = n[n[0].astype(str).str.startswith('RDBECOLITFC')]
E0 = n[[1, 4]].drop_duplicates().rename(columns={1: 'tf', 4: 'gene'})


def stats(targets, deg):
    d = deg.reindex(list(targets)).values
    return d.mean(), 100 * (d >= 3).mean()


def curveball(adj, trades):
    """adj: list of sets (targets of each TF); in-place curveball trades between random TF pairs."""
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


def analyse(E, fe, label):
    deg = E.groupby('gene').tf.nunique()
    tfs = sorted(E.tf.unique())
    targ = E.groupby('tf').gene.apply(set)
    fe = [t for t in fe if t in targ.index]
    T = set().union(*[targ[t] for t in fe])
    obs_mean, obs_pct = stats(T, deg)
    row = {'analysis': label, 'n_TFs': len(tfs), 'n_genes': len(deg), 'FE_TFs': ','.join(fe), 'n_FE_targets': len(T),
           'network_mean': deg.mean(), 'network_pct3': 100 * (deg >= 3).mean(),
           'obs_mean': obs_mean, 'obs_pct3': obs_pct}
    genes = deg.index.values
    # 1. uniform gene sets
    u = np.array([stats(rng.choice(genes, len(T), replace=False), deg) for _ in range(N_UNI)])
    # 2. degree-matched TF sets
    outdeg = targ.apply(len)
    others = outdeg.drop(fe)
    pools = [others.iloc[np.argsort(np.abs(others.values - outdeg[t]))[:6]].index.values for t in fe]
    m = []
    for _ in range(N_TF):
        pick = {rng.choice(p) for p in pools}
        m.append(stats(set().union(*[targ[t] for t in pick]), deg))
    m = np.array(m)
    # 3. degree-preserving network randomization (gene in-degree and TF out-degree fixed)
    adj = [set(targ[t]) for t in tfs]
    idx = [tfs.index(t) for t in fe]
    curveball(adj, 20000)  # burn-in
    r = []
    for _ in range(N_NET):
        curveball(adj, TRADES)
        r.append(stats(set().union(*[adj[i] for i in idx]), deg))
    r = np.array(r)
    for name, a in [('uniform', u), ('matchedTF', m), ('degreepres', r)]:
        row.update({f'{name}_mean': a[:, 0].mean(), f'{name}_mean_95hi': np.quantile(a[:, 0], 0.95),
                    f'{name}_pct3': a[:, 1].mean(),
                    f'{name}_p_mean': (1 + (a[:, 0] >= obs_mean).sum()) / (len(a) + 1),
                    f'{name}_p_pct3': (1 + (a[:, 1] >= obs_pct).sum()) / (len(a) + 1)})
    print(f"{label:28s} FE targets {len(T):4d}  obs mean {obs_mean:.2f} ({obs_pct:.1f}% >=3) | uniform {u[:,0].mean():.2f} "
          f"p={row['uniform_p_mean']:.4f} | matched TFs {m[:,0].mean():.2f} p={row['matchedTF_p_mean']:.4f} | "
          f"degree-preserving {r[:,0].mean():.2f} p={row['degreepres_p_mean']:.4f}")
    return row


rows = [analyse(E0, FE, 'all')]
for t in FE:
    rows.append(analyse(E0, [x for x in FE if x != t], f'FE set without {t}'))
for g in GLOBAL:
    rows.append(analyse(E0[E0.tf != g], FE, f'network without {g}'))
rows.append(analyse(E0[~E0.tf.isin(GLOBAL)], FE, 'network without 7 global TFs'))
pd.DataFrame(rows).to_csv(f'{OUT}/a31_tf_redundancy_null.csv', index=False)
print('missing global TF names:', [g for g in GLOBAL if g not in set(E0.tf)])
