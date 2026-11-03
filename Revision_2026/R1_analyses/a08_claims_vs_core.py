"""Which Results statements rely on fitness-enhancing genes outside the significant core?
For every gene the manuscript names as fitness-enhancing (Results l. 437-534 and Tables S17-S18), report
n conditions (of 46) with fitness > 0.5, and with fitness > 0.5 and t > 4, and whether it is a candidate
(>= 10) or core (>= 10, with t > 4 in >= 1 of them). Also re-computes Table 2 (regulatory redundancy) with core TFs only."""
import os
import numpy as np, pandas as pd
from common import *

ROOT = os.path.dirname(HERE)
SUPP = os.path.join(ROOT, 'external', 'Supp_Tables_1-19_Hidalgo_et_al.xlsx')
s19 = pd.read_excel(os.path.join(ROOT, 'Table_S19_revised.xlsx'), header=1)
s19 = s19[s19.bnumber.astype(str).str.match(r'^b\d{4}$')]
eco = pd.read_csv(os.path.join(REPO, 'Input_Data', 'All-genes-of-E.-coli-K-12-substr-MG1655.txt'), sep='\t')
alias = {}
for _, r in eco.dropna(subset=['Accession-1']).iterrows():
    alias[str(r['Gene Name']).lower()] = r['Accession-1']
for b, g in zip(s19.bnumber, s19.gene_name):
    alias.setdefault(str(g).lower(), b)
S = s19.set_index('bnumber')
ncol, tcol = 'sum_of_enhancing', 'sum_of_enhancing_significant (t > 4)'

claims = {
    'l.457 sensing/regulation genes': ['rpoS', 'fis', 'hfq', 'arcB', 'greA', 'cspC', 'rcsB'],
    'l.458 small regulatory RNAs': ['oxyS', 'arcZ', 'gcvB', 'dsrA'],
    'l.459 rsd': ['rsd'],
    'l.462 flagella': ['fliT', 'flgN'],
    'l.462 hibernation factors': ['rmf', 'ettA'],
    'l.463 ribosomal proteins': ['rpsG', 'rpsR'],
    'l.463 tRNAs (as stated)': ['selC', 'selU', 'valZ', 'thrT'],
    'l.507 TFs used for Table 2 / Fig. 6C': ['rbsR', 'gcvA', 'fis', 'slyA', 'rcsB'],
}
x17 = pd.read_excel(SUPP, sheet_name='Table_S17', header=None).dropna(how='all').iloc[1:]
x18 = pd.read_excel(SUPP, sheet_name='Table_S18', header=None).dropna(how='all').iloc[1:]
def genes_of(t):
    out = []
    for txt in t[1].astype(str):
        m = pd.Series(txt).str.extract(r'gene\s*:?\s*([A-Za-z0-9]+)')[0].iloc[0]
        out.append(m)
    return out
claims['Table S17 high-confidence predictions'] = genes_of(x17)
claims['Table S18 low-confidence / no prediction'] = genes_of(x18)

rows = []
for claim, genes in claims.items():
    for g in genes:
        b = alias.get(str(g).lower())
        if b in S.index:
            n, nt = int(S.at[b, ncol]), int(S.at[b, tcol])
            status = 'core' if (n >= 10 and nt >= 1) else ('candidate only' if n >= 10 else 'not a candidate')
        else:
            n = nt = np.nan; status = 'no fitness data'
        rows.append({'claim': claim, 'gene': g, 'bnumber': b, 'n_cond_fitness>0.5': n, 'n_cond_t>4': nt, 'status': status})
res = pd.DataFrame(rows)
res.to_csv(f'{OUT}/a08_claims_vs_core.csv', index=False)
pd.set_option('display.width', 200)
print(res.to_string(index=False))
print('\nSummary by claim:')
print(res.groupby('claim').status.value_counts().unstack(fill_value=0).to_string())

# rpoS: "fitness-enhancing in seven out of ten growth conditions that do not involve stress" (shared conditions)
fit, _ = load_fitness()
tt = pd.read_csv(os.path.join(ROOT, 'Input_Data_Price2018', 'fit_t.tab'), sep='\t').set_index('sysName')
ft = pd.read_csv(os.path.join(REPO, 'Input_Data', 'fit_organism_Keio.tsv'), sep='\t').set_index('sysName')
b = 'b2741'
r = []
for c in CONDS:
    col = FIT_COL[c]
    pair = ['set2IT045 LB', 'set2IT069 LB'] if c == 'LB' else [col, ft.columns[list(ft.columns).index(col) + 1]]
    t = tt.loc[b, pair].sum() / np.sqrt(3)
    r.append({'condition': c, 'fitness': round(fit.at[b, c], 2), 't': round(t, 1)})
r = pd.DataFrame(r)
print('\nrpoS in the 11 shared conditions:'); print(r.to_string(index=False))
me = r[r.condition != 'LB']
print(f'rpoS fitness > 0.5 in {(me.fitness > 0.5).sum()}/10 minimal-medium conditions; '
      f'with t > 4: {((me.fitness > 0.5) & (me.t > 4)).sum()}/10')

# ---- Table 2 with the five TFs vs the core TFs ---------------------------------------------------
net = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'NetworkRegulatorGene.tsv'), sep='\t',
                  comment='#', header=None).iloc[:, [1, 4, 5]]
net.columns = ['TF', 'gene', 'effect']
net = net.dropna(subset=['TF', 'gene'])
def table2(tfs):
    targets = set(net.loc[net.TF.str.lower().isin([t.lower() for t in tfs]), 'gene'])
    def stats_for(genes):
        sub = net[net.gene.isin(genes)]
        ntf = sub.groupby('gene').TF.nunique()
        red = sub.groupby(['gene', 'effect']).TF.nunique().groupby('gene').max() >= 2
        return {'n_genes': len(ntf), '%_1_TF': 100 * (ntf == 1).mean(), '%_2_TFs': 100 * (ntf == 2).mean(),
                '%_>=3_TFs': 100 * (ntf >= 3).mean(), '%_redundant (same effect by >=2 TFs)': 100 * red.reindex(ntf.index, fill_value=False).mean(),
                'mean_TFs_per_gene': ntf.mean()}
    return stats_for(targets)
t2 = pd.DataFrame({'TRN (all regulated genes)': table2(net.TF.unique()),
                   'five TFs (RbsR, GcvA, Fis, SlyA, RcsB)': table2(['RbsR', 'GcvA', 'Fis', 'SlyA', 'RcsB']),
                   'core TFs (RbsR, GcvA, Fis)': table2(['RbsR', 'GcvA', 'Fis'])}).T
t2.round(1).to_csv(f'{OUT}/a08_table2_core_TFs.csv')
print('\nTable 2 recomputed (RegulonDB network in the repository):'); print(t2.round(1).to_string())
