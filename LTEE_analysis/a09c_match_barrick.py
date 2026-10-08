"""Match every ALEdb "LTEE ARA" clone to its curated breseq GenomeDiff in barricklab/LTEE-Ecoli (LTEE-clone-curated),
to identify the populations (A1-A6) and their mutator status, and resolve what every DEL with within= sits in.

GenomeDiffs: LTEE-clone-curated/ (copied from barricklab/LTEE-Ecoli, CC0; set LTEE_GD for another copy).
Matching: Jaccard index of (mutation type, position) sets; the best-matching GenomeDiff file per clone.
"""
import os, glob, collections
import pandas as pd
from common import *

G = os.environ.get('LTEE_GD', os.path.join(HERE, 'LTEE-clone-curated'))
MUT = {'SNP', 'SUB', 'DEL', 'INS', 'MOB', 'AMP', 'CON', 'INV'}
gd, meta, recs_all = {}, {}, {}
for f in glob.glob(os.path.join(G, '*.gd')):
    s, m, recs = set(), {}, {}
    for line in open(f, encoding='utf-8'):
        if line.startswith('#='):
            k, _, v = line[2:].strip().partition('\t'); m[k] = v
        else:
            p = line.rstrip('\n').split('\t')
            if p[0] in MUT:
                s.add((p[0], int(p[4]))); recs[p[1]] = p
    n = os.path.basename(f)[:-3]
    gd[n], meta[n], recs_all[n] = s, m, recs

d = pd.read_csv(os.path.join(INPUT, 'Proj_LTEE_Exp_LTEE_ARA_ExpID350_mut.csv'))
d['pos'] = d.Position.astype(str).str.replace(',', '').str.extract(r'(\d+)')[0].astype(int)
rows = []
for c in [c for c in d.columns if c.startswith('LTEE')]:
    s = set(zip(d.loc[d[c].notna(), 'Mutation Type'], d.loc[d[c].notna(), 'pos']))
    best = max(gd, key=lambda n: len(s & gd[n]) / len(s | gd[n]))
    rows.append({'aledb_clone': c.replace('LTEE ARA ', '').replace(' R1', ''), 'n_aledb': len(s), 'barrick_gd': best,
                 'n_gd': len(gd[best]), 'shared': len(s & gd[best]), 'jaccard': round(len(s & gd[best]) / len(s | gd[best]), 3),
                 'population': meta[best].get('POPULATION'), 'mutator_status': meta[best].get('MUTATOR_STATUS')})
r = pd.DataFrame(rows)
r.to_csv(f'{OUT}/a09c_aledb_to_barrick_clones.csv', index=False)
r['A'] = r.aledb_clone.str[:2]
print(r.groupby('A').agg(population=('population', lambda x: ','.join(sorted(set(x)))),
                         mutator=('mutator_status', lambda x: ','.join(sorted(set(x)))),
                         clones=('aledb_clone', 'size'), jaccard_median=('jaccard', 'median'),
                         aledb_only_mutations=('n_aledb', lambda x: int((x - r.loc[x.index, 'shared']).sum()))).to_string())

mm = pd.DataFrame([{'population': m.get('POPULATION'), 'generation': int(m.get('TIME', 0)), 'status': m.get('MUTATOR_STATUS')}
                   for m in meta.values()])
st = mm.groupby(['population', 'status']).generation.agg(['min', 'max', 'count']).reset_index()
st.to_csv(f'{OUT}/a09c_barrick_mutator_status.csv', index=False)
print('\nMutator status of all curated LTEE clones (Barrick):\n', st.to_string(index=False))

res = collections.defaultdict(set)
for n, recs in recs_all.items():
    if not n.startswith('Ara+'):
        continue
    for p in recs.values():
        w = [a for a in p[6:] if a.startswith('within=')]
        if p[0] == 'DEL' and w:
            t = recs.get(w[0].split('=')[1].split(':')[0])
            res[(int(p[4]), int(p[5]))].add(f"{t[0]} {t[4]} {t[5]}" + (f" x{t[6]}" if t[0] == 'AMP' else '') if t else 'missing')
w = pd.DataFrame([{'position': k[0], 'size_bp': k[1], 'within': '; '.join(sorted(v)),
                   'inside_amplification': any(x.startswith('AMP') for x in v)} for k, v in sorted(res.items())])
w.to_csv(f'{OUT}/a09c_deletions_within.csv', index=False)
print('\nDeletions with within= in Ara+ GenomeDiffs:\n', w.to_string(index=False))
