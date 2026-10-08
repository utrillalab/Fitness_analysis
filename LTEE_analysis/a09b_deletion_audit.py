"""Audit of large (>=5 kb) deletions in the ALEdb LTEE file: which are novel junctions of amplifications?

For every deletion: amplifications with the same start coordinate carried by the same clone, other mutations lying
inside the "deleted" span in the same clone (impossible if the span were really lost), genes in the span with no
RB-TnSeq data (likely essential) and genes with glucose fitness <= -2.
Also traces the source event of each of the 11 fitness-enhancing LoF genes of the manuscript (original fitness).
"""
import os, re
import pandas as pd
from common import *

d = pd.read_csv(os.path.join(INPUT, 'Proj_LTEE_Exp_LTEE_ARA_ExpID350_mut.csv'))
cl = [c for c in d.columns if c.startswith('LTEE')]
d['pos'] = d.Position.astype(str).str.replace(',', '').str.extract(r'(\d+)')[0].astype(float)


def bp(ch):
    m = re.match(r'^Δ?([\d,]+) bp', str(ch).replace(',', ''))
    return int(m.group(1)) if m else None


d['size'] = d['Sequence Change'].apply(bp)
d.loc[d['Mutation Type'] == 'SUB', 'size'] = d.loc[d['Mutation Type'] == 'SUB', 'Sequence Change'].str.extract(r'^([\d,]+) bp')[0].str.replace(',', '').astype(float)
d['cset'] = d[cl].notna().apply(lambda r: {c for c, v in zip(cl, r) if v}, axis=1)
d['genes'] = d['Gene (Scrollable)'].fillna('').apply(lambda s: [g.strip(' []') for g in s.split(',') if g.strip(' []')])
fit, names = load_fitness(); glc = fit.Glucose
eco = pd.read_csv(os.path.join(REPO, 'Input_Data', 'All-genes-of-E.-coli-K-12-substr-MG1655.txt'), sep='\t')
n2b = dict(zip(eco['Gene Name'].astype(str), eco['Accession-1'])); n2b.update({str(n): b for b, n in names.items()})
amp = d[d['Mutation Type'] == 'AMP']
short = lambda cs: ', '.join(sorted(c.replace('LTEE ARA ', '').replace(' R1', '') for c in cs))

rows = []
for i, r in d[d['Mutation Type'].isin(['DEL', 'SUB']) & (d['size'].fillna(0) >= 5000)].iterrows():
    s, e = r.pos, r.pos + r['size']
    same = [a for _, a in amp.iterrows() if abs(a.pos - s) <= 10 and (r.cset & a.cset)]
    ins = d[(d.pos > s + 50) & (d.pos < e - 50) & (d.index != i) & (d['Mutation Type'] != 'AMP')]
    inside = [f"{x['Gene (Scrollable)']} {x['Mutation Type']}" for _, x in ins.iterrows() if x.cset & r.cset]
    gb = [(g, n2b.get(g)) for g in r.genes]
    rows.append({'position': int(s), 'type': r['Mutation Type'], 'size_bp': int(r['size']), 'n_genes': len(r.genes),
                 'clones': short(r.cset), 'contains_rpoS': 'rpoS' in r.genes,
                 'same-start AMP in same clone': '; '.join(f"{a['Sequence Change']} ({len(r.cset & a.cset)}/{len(r.cset)} clones)" for a in same),
                 'mutations inside span, same clone': '; '.join(inside[:5]),
                 'genes without RB-TnSeq data (likely essential)': ', '.join(g for g, b in gb if b and pd.isna(glc.get(b)) and not g.startswith('ins')),
                 'genes with glucose fitness <= -2': ', '.join(g for g, b in gb if b in glc.index and glc[b] <= -2),
                 'verdict': 'amplification junction' if same else 'deletion'})
audit = pd.DataFrame(rows).sort_values('position')
audit.to_csv(f'{OUT}/a09b_large_deletions_audit.csv', index=False)
print(audit[audit.verdict == 'amplification junction'].drop(columns='verdict').to_string(index=False))

# source of the manuscript's 11 enhancing LoF genes (original repo fitness, all deletions counted)
orig = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'genes_categories', 'genesfitness_Glucose.csv')).set_index('bnum')
genes = pd.read_csv(f'{OUT}/a09_ltee_lof_genes.csv')
src = []
for b, f in orig[orig.Categ.str.startswith('Enh')].fitval.items():
    g = names.get(b)
    ev = d[d.genes.apply(lambda gs: g in gs) & d['Mutation Type'].isin(['DEL', 'SUB', 'MOB', 'INS', 'SNP'])]
    if not len(ev):
        continue
    src.append({'gene': g, 'bnumber': b, 'fitness_original': f, 'fitness_corrected': glc.get(b),
                'events': '; '.join(f"{x['Mutation Type']} {x['Sequence Change']} {str(x.Details)[:20]}" for _, x in ev.iterrows()),
                'real LoF after audit': b in set(genes.bnumber)})
src = pd.DataFrame(src)
src.to_csv(f'{OUT}/a09b_manuscript_enhancing_sources.csv', index=False)
print('\nSources of enhancing genes hit in the LTEE (original fitness):\n', src.round(2).to_string(index=False))
