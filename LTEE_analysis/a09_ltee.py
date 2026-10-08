"""R1-M6 / R2-M3: loss-of-function (LoF) mutations in the LTEE vs fitness categories in glucose.

Data: ALEdb (Phaneuf et al. 2019), experiment "LTEE ARA" (ExpID 350; Tenaillon et al. 2016), exported to
`Input_Data/Proj_LTEE_Exp_LTEE_ARA_ExpID350_mut.csv`. 679 mutations in 106 clones from 6 populations
(A1-A6). Each clone column is non-empty when the clone carries the mutation.

LoF (Methods of the manuscript: frameshift, nonsense, large deletion, IS disruption):
  nonsense SNP; coding indel whose length is not a multiple of 3; IS insertion inside a coding region;
  deletion (or SUB) that removes one or more whole/partial genes (IS-mediated or large).
Not LoF: missense, synonymous, in-frame indels, intergenic, pseudogene, amplifications, inversion, conversion,
and "deletions" starting where an amplification carried by the same clone starts (novel junction of an amplification in the same clone).

Fitness: glucose, mean of IT003/IT004 (corrected; common.load_fitness). t from Price et al. 2018 (Input_Data/fit_t_glucose_Price2018.tsv: fit_t.tab, glucose columns).
Universe: genes with glucose fitness that are in the EcoCyc K-12 gene table (needed for gene length/position).
REL606 gene names are matched to K-12 b-numbers by name; IS elements (ins*) and ECB_ loci are not mapped.

Tests (gene level, union of genes hit by >=1 LoF event):
  1. hypergeometric per category vs the universe (the test R2 asks for);
  2. event-level permutation (10,000): every single-gene LoF event is replaced by a random gene drawn with
     probability proportional to its length (mutational target size) and every multi-gene deletion by a random
     block of the same number of consecutive K-12 genes (deletions remove contiguous blocks, not independent genes);
  3. Mann-Whitney of glucose fitness, LoF genes vs the rest of the universe.
Subsets: all LoF; point-like LoF only (no multi-gene deletions); parallel (LoF in >=2 populations);
persistent (same mutation seen at >=2 time points). Also missense/synonymous genes as comparators.
Hypermutators: A1-A6 = Ara+1..Ara+6 (matched to barricklab/LTEE-Ecoli GenomeDiffs); ALEdb has no point-mutator
clones, but Ara+1 clones >= 5,000 gen are IS-mutators -> sensitivity subset without them.
"""
import os, re, sys, glob
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from common import *

rng = np.random.default_rng(1)
NPERM = 10000
# Which populations: 'plus' (ALEdb "LTEE ARA", Ara+1..6; the set used in the manuscript), 'minus' (ALEdb "LTEE",
# Ara-1..6; independent replicate) or 'all'. Outputs: a09_* for plus, a09_minus_* and a09_all_* otherwise.
SET = sys.argv[1] if len(sys.argv) > 1 else 'plus'
SUF = '' if SET == 'plus' else f'_{SET}'
ALE_FILES = {'plus': ('Ara+', os.path.join(INPUT, 'Proj_LTEE_Exp_LTEE_ARA_ExpID350_mut.csv')),
             'minus': ('Ara-', os.path.join(INPUT, 'ALEdb_LTEE_AraMinus_mut.csv'))}

# ---------------------------------------------------------------------------------------------------
# 1. fitness, t, categories, gene table
fit, names = load_fitness()
glc = fit['Glucose']
tt = pd.read_csv(os.path.join(INPUT, 'fit_t_glucose_Price2018.tsv'), sep='\t').set_index('sysName')
t_glc = (tt[['set1IT003 D-Glucose (C)', 'set1IT004 D-Glucose (C)']].sum(axis=1) / np.sqrt(2) / np.sqrt(1.5)).reindex(glc.index)
cat = pd.Series(classify(glc), index=glc.index)
cat_t = cat.where((cat == 'mean-effect') | (cat == 'NA') | (t_glc.abs() > 4), 'mean-effect')

eco = pd.read_csv(os.path.join(REPO, 'Input_Data', 'All-genes-of-E.-coli-K-12-substr-MG1655.txt'), sep='\t')
eco = eco.dropna(subset=['Accession-1', 'Left-End-Position', 'Right-End-Position'])
eco = eco[eco['Accession-1'].str.match(r'^b\d{4}$')].drop_duplicates('Accession-1')
eco['len'] = (eco['Right-End-Position'] - eco['Left-End-Position']).abs() + 1
eco = eco.sort_values('Left-End-Position').reset_index(drop=True)
universe = [b for b in eco['Accession-1'] if b in glc.index and pd.notna(glc[b])]
U = set(universe)

# 108 recurrent fitness-enhancing candidates (fitness > 0.5 in >= 10 of 46 conditions); core = also t > 4 in >= 1 (71 genes)
core = pd.read_csv(os.path.join(INPUT, 'enhancing_candidates_108_core.csv')).set_index('bnumber')
net = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'NetworkRegulatorGene.tsv'), sep='\t', comment='#', header=None)
reg_lower = {str(r).lower() for r in net[1].dropna()}
product = eco.set_index('Accession-1').Product.astype(str)

name2b = {}
for b, n in eco[['Accession-1', 'Gene Name']].values:
    name2b.setdefault(str(n), b)
for b, n in names.items():
    if pd.notna(n):
        name2b[str(n)] = b          # fitness-file names take precedence

# ---------------------------------------------------------------------------------------------------
# 2. parse mutations
parts = []
for key in (['plus', 'minus'] if SET == 'all' else [SET]):
    prefix, path = ALE_FILES[key]
    x = pd.read_csv(path)
    x.columns = [re.sub('<.*', '', c).strip() for c in x.columns]
    # clone columns: 'LTEE ARA A1 F500 I1 R1' or 'A1 F500 I1 R1' -> 'Ara+1 F500 I1'
    ren = {c: re.sub(r'^(?:LTEE ARA )?A(\d+) (F\d+ I\d+).*$', lambda m: f'{prefix}{m.group(1)} {m.group(2)}', c)
           for c in x.columns if re.match(r'^(?:LTEE ARA )?A[1-9]\d* F\d+ I\d+', c)}
    x = x.rename(columns=ren)[['Reference Seq', 'Position', 'Mutation Type', 'Sequence Change', 'Gene (Scrollable)', 'Details'] + list(ren.values())]
    parts.append(x)
d = pd.concat(parts, ignore_index=True)
clones = [c for c in d.columns if c.startswith('Ara')]
d = d[d[clones].notna().any(axis=1)].reset_index(drop=True)
cinfo = pd.DataFrame({'clone': clones})
cinfo[['population', 'gen']] = cinfo.clone.str.extract(r'^(Ara[+-]\d) F(\d+) ')
# A1-A6 of the "LTEE ARA" export are Ara+1..6 and those of the "LTEE" export are Ara-1..6: every ALEdb clone matches
# its curated GenomeDiff in barricklab/LTEE-Ecoli (a09c_match_barrick.py). ALEdb keeps only clones that are not
# point-mutators; some are IS-mutators (Ara+1 from 5,000 gen; one Ara-5 clone at 30,000), see section 2b.
cinfo['gen'] = cinfo.gen.astype(int)


def indel_len(change):
    s = str(change).replace(',', '')
    m = re.match(r'^\((\w+)\)(\d+)→(\d+)$', s)
    if m:
        return len(m.group(1)) * (int(m.group(3)) - int(m.group(2)))
    m = re.match(r'^\+([ACGT]+)$', s)
    if m:
        return len(m.group(1))
    m = re.match(r'^\+(\d+) bp', s)
    if m:
        return int(m.group(1))
    m = re.match(r'^Δ(\d+) bp', s)
    if m:
        return -int(m.group(1))
    return None


def mut_class(r):
    typ, det, ch = r['Mutation Type'], str(r['Details']), str(r['Sequence Change'])
    coding = det.startswith('coding')
    if typ == 'SNP':
        m = re.match(r'^([A-Z\*])\d+([A-Z\*])', det)
        if not m:
            return 'other'
        if m.group(2) == '*' and m.group(1) != '*':
            return 'nonsense'
        return 'synonymous' if m.group(1) == m.group(2) else 'missense'
    if typ in ('INS', 'DEL') and coding:
        n = indel_len(ch)
        if n is None:
            return 'other'
        return 'frameshift' if n % 3 else 'inframe_indel'
    if typ == 'DEL' and not det.startswith(('intergenic', 'pseudogene', 'noncoding')):  # incl. 'between IS' deletions
        return 'deletion'
    if typ == 'SUB':
        m = re.match(r'^([\d,]+) bp→', ch)
        if m and int(m.group(1).replace(',', '')) >= 50 and not det.startswith('intergenic'):
            return 'deletion'
        return 'other'
    if typ == 'MOB':
        return 'IS_insertion' if coding else 'other'
    return 'other'


LOF = {'nonsense', 'frameshift', 'deletion', 'IS_insertion'}
d['class'] = d.apply(mut_class, axis=1)
d['genes'] = d['Gene (Scrollable)'].fillna('').apply(lambda s: [g.strip(' []') for g in s.split(',') if g.strip(' []')])
d['genes'] = d['genes'].apply(lambda gs: list(dict.fromkeys(gs)))

# "Deletions" that start at the same coordinate as an amplification carried by the SAME clone are annotation
# artefacts of that amplification, not real deletions: a clone cannot carry 2-3 copies and 0 copies of the same
# span; in 3 of 6 cases other mutations of the same clone lie inside the "deleted" span; and spans contain essential
# or auxotrophy genes (the 60,769-bp event with rpoS would remove alaS, ispD, ispF, ftsB and cysCND).
# See a09b_deletion_audit.py.
def bp(ch):
    m = re.match(r'^Δ?([\d,]+) bp', str(ch))
    return int(m.group(1).replace(',', '')) if m else None


d['pos'] = d.Position.astype(str).str.replace(',', '').str.extract(r'(\d+)')[0].astype(float)
d['cset'] = d[clones].notna().apply(lambda row: {c for c, v in zip(clones, row) if v}, axis=1)
amps = d[d['Mutation Type'] == 'AMP']
for i, r in d[d['class'] == 'deletion'].iterrows():
    if any(abs(r.pos - a.pos) <= 10 and (r.cset & a.cset) for _, a in amps.iterrows()):
        d.at[i, 'class'] = 'amplification_artefact'

# Cross-check with the curated breseq GenomeDiff files (barricklab/LTEE-Ecoli, LTEE-clone-curated), where these
# deletions carry within=<AMP id>:<copy>, i.e. they happen inside one copy of the amplification. ALEdb drops that
# attribute.
# Curated breseq GenomeDiffs of barricklab/LTEE-Ecoli (CC0), copied to LTEE-clone-curated/ (commit 2307710, 2026-06-23).
GD_DIR = os.environ.get('LTEE_GD', os.path.join(HERE, 'LTEE-clone-curated'))
gd_dir = GD_DIR
if os.path.isdir(gd_dir):
    within_amp = set()
    for f in os.listdir(gd_dir):
        if f.endswith('.gd') and f.split('_')[0] in set(cinfo['population']):
            recs = {}
            for line in open(os.path.join(gd_dir, f), encoding='utf-8'):
                p = line.rstrip('\n').split('\t')
                if p[0] in ('SNP', 'SUB', 'DEL', 'INS', 'MOB', 'AMP', 'CON', 'INV'):
                    recs[p[1]] = p
            for p in recs.values():
                w = [a.split('=')[1].split(':')[0] for a in p[6:] if a.startswith('within=')]
                if p[0] == 'DEL' and w and recs.get(w[0], [''])[0] == 'AMP':
                    within_amp.add((f.split('_')[0], int(p[4])))
    # The GenomeDiffs are the primary criterion: the ALEdb Ara- export has no AMP rows at all, so the same-start
    # rule above cannot see, e.g., the 26,603-bp "deletion" (fhlA..mutS) of Ara-5, which is within an AMP copy.
    pop_of = dict(zip(cinfo.clone, cinfo['population']))
    rule = set(d.index[d['class'] == 'amplification_artefact'])
    for i, r in d[d['class'] == 'deletion'].iterrows():
        if any((pop_of[c], int(r.pos)) in within_amp for c in r.cset):
            d.at[i, 'class'] = 'amplification_artefact'
    gdset = set(d.index[d['class'] == 'amplification_artefact'])
    print(f'GenomeDiff check: deletions inside an AMP copy = {sorted(d.loc[sorted(gdset), "pos"].astype(int))}; '
          f'same-start rule found {len(rule)} of {len(gdset)} rows'
          f'{"" if rule <= gdset else "; rule flagged rows not confirmed by GenomeDiff: " + str(sorted(rule - gdset))}')
    # 2b. mutator status of each ALEdb clone = that of its best-matching curated GenomeDiff
    gsets, gstat = {}, {}
    for f in glob.glob(os.path.join(gd_dir, 'Ara*.gd')):
        st, ms = None, set()
        for line in open(f, encoding='utf-8'):
            if line.startswith('#=MUTATOR_STATUS'):
                st = line.strip().split('	')[1]
            pp = line.split('	')
            if pp[0] in ('SNP', 'SUB', 'DEL', 'INS', 'MOB', 'AMP', 'CON', 'INV'):
                ms.add((pp[0], int(pp[4])))
        gsets[f] = ms; gstat[f] = st
    status = {}
    for c in clones:
        sc = set(zip(d.loc[d[c].notna(), 'Mutation Type'], d.loc[d[c].notna(), 'pos'].astype(int)))
        best = max(gsets, key=lambda f: len(sc & gsets[f]) / max(len(sc | gsets[f]), 1))
        status[c] = gstat[best]
    cinfo['mutator_status'] = cinfo.clone.map(status)
else:
    print('GenomeDiff check skipped (LTEE-Ecoli repo not found); IS-mutator clones set from Barrick metadata by hand')
    cinfo['mutator_status'] = np.where((cinfo['population'] == 'Ara+1') & (cinfo.gen >= 5000), 'IS-mutator', 'non-mutator')
print('Mutator status of the ALEdb clones:', cinfo.groupby(['population', 'mutator_status']).size().to_dict())
d['pops'] = d[clones].notna().apply(lambda row: sorted({cinfo['population'][i] for i, v in enumerate(row) if v}), axis=1)
d['n_timepoints'] = d[clones].notna().apply(lambda row: len({(cinfo['population'][i], cinfo.gen[i]) for i, v in enumerate(row) if v}), axis=1)
d['bnums'] = d.genes.apply(lambda gs: [name2b[g] for g in gs if g in name2b])
d['block'] = (d['class'] == 'deletion') & (d.genes.str.len() > 1)

print('Mutation classes:\n', d['class'].value_counts().to_string())
lof = d[d['class'].isin(LOF)].copy()
allg = sorted({g for gs in lof.genes for g in gs})
mapped = sorted({g for g in allg if g in name2b})
unm = [g for g in allg if g not in name2b]
print(f'\nLoF events: {len(lof)}; gene names hit: {len(allg)}; mapped to b-number: {len(mapped)}; '
      f'unmapped: {len(unm)} (IS elements {sum(g.startswith("ins") for g in unm)}, ECB_ loci '
      f'{sum(g.startswith("ECB_") for g in unm)}, other: {[g for g in unm if not g.startswith(("ins", "ECB_"))]})')

# ---------------------------------------------------------------------------------------------------
# 3. hypermutator check
cnt = d[clones].notna().sum()
cinfo['n_mut'] = cnt.values
cinfo['n_snp'] = d.loc[d['Mutation Type'] == 'SNP', clones].notna().sum().values
fitrate = stats.linregress(cinfo.gen, cinfo.n_mut)
cinfo['resid'] = cinfo.n_mut - (fitrate.intercept + fitrate.slope * cinfo.gen)
repair = r'\b(mutS|mutL|mutH|mutT|mutY|mutM|uvrA|uvrB|uvrD|dnaQ|polB|recA)\b'
rep = d[d['Gene (Scrollable)'].fillna('').str.contains(repair)]
popsum = []
for p, g in cinfo.groupby('population'):
    s = d[d[g.clone.tolist()].notna().any(axis=1)]
    popsum.append({'population': p, 'clones': len(g), 'generations': f'{g.gen.min()}-{g.gen.max()}',
                   'max_mutations_per_clone': int(g.n_mut.max()), 'max_SNPs_per_clone': int(g.n_snp.max()),
                   'mutations_per_1000_gen': round(np.polyfit(g.gen, g.n_mut, 1)[0] * 1000, 2),
                   'IS150_%': round(100 * s['Sequence Change'].fillna('').str.startswith('IS150').mean(), 1),
                   'DNA_repair_genes_hit': '; '.join(f"{r['Gene (Scrollable)'][:25]} {r['Mutation Type']} {r['Details']}"
                                                     for _, r in rep.iterrows() if r[g.clone.tolist()].notna().any())})
popsum = pd.DataFrame(popsum)
popsum.to_csv(f'{OUT}/a09{SUF}_populations_mutator_check.csv', index=False)
cinfo.to_csv(f'{OUT}/a09{SUF}_clones.csv', index=False)
print('\nPopulations:\n', popsum.to_string(index=False))

# ---------------------------------------------------------------------------------------------------
# 4. tests
gpos = {b: i for i, b in enumerate(eco['Accession-1'])}
all_b = eco['Accession-1'].values
w = eco['len'].values / eco['len'].sum()


def genes_of(events):
    return sorted({b for bs in events.bnums for b in bs} & U)


def null_sets(events, n=NPERM):
    """Event-level null: same number of single-gene events (length-weighted) and same-size deletion blocks."""
    singles = int((~events.block).sum())
    sizes = events.loc[events.block, 'genes'].str.len().values
    out = []
    for _ in range(n):
        s = set(all_b[rng.choice(len(all_b), size=singles, p=w)])
        for k in sizes:
            i = rng.integers(0, len(all_b) - k)
            s.update(all_b[i:i + k])
        out.append(list(s & U))
    return out


def summarize(label, events, catser, null=None):
    genes = genes_of(events)
    c = catser.reindex(universe)
    N = len(universe); n = len(genes)
    row = {'set': label, 'events': len(events), 'genes_in_universe': n}
    for k in CATS:
        K = int((c == k).sum()); x = int((c.reindex(genes) == k).sum())
        row[f'{k}_n'] = x
        row[f'{k}_%'] = 100 * x / max(n, 1)
        row[f'{k}_background_%'] = 100 * K / N
        row[f'{k}_fold'] = (x / max(n, 1)) / (K / N)
        row[f'{k}_p_hyper_enriched'] = stats.hypergeom.sf(x - 1, N, K, n)
        row[f'{k}_p_hyper_depleted'] = stats.hypergeom.cdf(x, N, K, n)
        if null is not None:
            nf = np.array([(c.reindex(s) == k).mean() for s in null])
            row[f'{k}_null_%_median'] = 100 * np.median(nf)
            row[f'{k}_p_perm_enriched'] = (1 + (nf >= x / n).sum()) / (1 + len(nf))
            row[f'{k}_p_perm_depleted'] = (1 + (nf <= x / n).sum()) / (1 + len(nf))
    f_in = glc.reindex(genes); f_out = glc.reindex(sorted(U - set(genes)))
    row['median_fitness_LoF'] = f_in.median(); row['median_fitness_rest'] = f_out.median()
    row['p_MWU_fitness'] = stats.mannwhitneyu(f_in, f_out).pvalue
    if null is not None:
        nm = np.array([glc.reindex(s).mean() for s in null])
        row['mean_fitness_LoF'] = f_in.mean(); row['mean_fitness_null_median'] = np.median(nm)
        row['p_perm_mean_fitness_higher'] = (1 + (nm >= f_in.mean()).sum()) / (1 + len(nm))
    return row


is_mut = set(cinfo.loc[cinfo.mutator_status != 'non-mutator', 'clone'])
subsets = {
    'all LoF': lof,
    'point-like LoF (no multi-gene deletions)': lof[~lof.block],
    'LoF in >=2 populations (parallel)': None,
    'LoF seen at >=2 time points (persistent)': lof[lof.n_timepoints >= 2],
    'LoF without IS-mutator clones': lof[lof.cset.apply(lambda cs: bool(cs - is_mut))],
}
gene_pops = {}
for _, r in lof.iterrows():
    for b in r.bnums:
        gene_pops.setdefault(b, set()).update(r.pops)
par = {b for b, p in gene_pops.items() if len(p) >= 2}
subsets['LoF in >=2 populations (parallel)'] = lof[lof.bnums.apply(lambda bs: any(b in par for b in bs))].assign(
    bnums=lambda x: x.bnums.apply(lambda bs: [b for b in bs if b in par]))

rows = []
for label, ev in subsets.items():
    null = null_sets(ev) if label != 'LoF in >=2 populations (parallel)' else None
    for cl, cs in [('fitness only', cat), ('fitness + |t|>4', cat_t)]:
        r = summarize(label, ev, cs, null)
        r['categories'] = cl
        rows.append(r)
for cls in ['missense', 'synonymous']:
    ev = d[(d['class'] == cls)]
    for cl, cs in [('fitness only', cat), ('fitness + |t|>4', cat_t)]:
        r = summarize(f'comparator: {cls}', ev, cs); r['categories'] = cl; rows.append(r)
res = pd.DataFrame(rows)
res = res[['set', 'categories'] + [c for c in res.columns if c not in ('set', 'categories')]]
res.to_csv(f'{OUT}/a09{SUF}_ltee_lof_tests.csv', index=False)
show = ['set', 'categories', 'events', 'genes_in_universe', 'mean-effect_%', 'mean-effect_background_%',
        'mean-effect_p_hyper_enriched', 'mean-effect_null_%_median', 'mean-effect_p_perm_enriched',
        'enhancing_n', 'enhancing_fold', 'enhancing_p_hyper_enriched', 'enhancing_p_perm_enriched',
        'essential_n', 'essential_p_hyper_depleted', 'important_n', 'p_MWU_fitness', 'p_perm_mean_fitness_higher']
print('\n', res[show].round(4).to_string(index=False))

# ---------------------------------------------------------------------------------------------------
# 5. gene table and the genes named in the text
gt = []
for b in sorted({b for bs in lof.bnums for b in bs}):
    ev = lof[lof.bnums.apply(lambda bs: b in bs)]
    gt.append({'bnumber': b, 'gene': names.get(b, product.get(b)), 'product': product.get(b),
               'LoF_classes': '; '.join(sorted(set(ev['class']))), 'only_in_multigene_deletions': bool(ev.block.all()),
               'populations': ','.join(sorted(gene_pops[b])), 'n_populations': len(gene_pops[b]),
               'glucose_fitness': glc.get(b), 'glucose_t': t_glc.get(b), 'category': cat.get(b, 'NA'),
               'category_t4': cat_t.get(b, 'NA'), 'in_universe': b in U,
               'enhancing_candidate_108': b in core.index, 'enhancing_core': bool(core.significant_core.get(b, False)),
               'RegulonDB_regulator': str(names.get(b, '')).lower() in reg_lower,
               'transporter (EcoCyc product)': bool(re.search(r'transport|permease|porin|symporter|antiporter|PTS',
                                                               str(product.get(b)), re.I))})
gt = pd.DataFrame(gt).sort_values(['category', 'glucose_fitness'], ascending=[True, False])
gt.to_csv(f'{OUT}/a09{SUF}_ltee_lof_genes.csv', index=False)
u = gt[gt.in_universe]
print(f'\nLoF genes mapped: {len(gt)}, with glucose fitness: {len(u)}, without: {len(gt) - len(u)}')
print(u.category.value_counts().to_string())
print('\nEnhancing LoF genes:\n', u[u.category == 'enhancing'][['gene', 'glucose_fitness', 'glucose_t', 'category_t4', 'LoF_classes',
      'n_populations', 'enhancing_candidate_108', 'enhancing_core']].round(2).to_string(index=False))
for g in ['rpoS', 'yebK', 'allR', 'malT', 'rbsR', 'nadR', 'pykF', 'spoT']:
    b = name2b.get(g)
    print(f'  {g} ({b}): fitness {glc.get(b, np.nan):.2f}, t {t_glc.get(b, np.nan):.1f}, LoF in LTEE: {b in set(gt.bnumber)}'
          f'{", pops " + ",".join(sorted(gene_pops[b])) if b in gene_pops else ""}')
for lbl, col in [('regulators', 'RegulonDB_regulator'), ('transporters', 'transporter (EcoCyc product)')]:
    s = u[u[col]]
    print(f'\n{lbl}: {len(s)}; categories: {s.category.value_counts().to_dict()}')

# ---------------------------------------------------------------------------------------------------
# 6. try to reproduce the manuscript counts (459 genes: 423 mean-effect, 13 essential, 12 important, 11 enhancing)
orig = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'genes_categories', 'genesfitness_Glucose.csv')).set_index('bnum')
# optional: published Supplementary Table S1 (not in the repo); set TABLE_S1 to its path to include it
S1 = os.environ.get('TABLE_S1', '')
s1 = (pd.read_excel(S1, sheet_name='Table S1', header=None, skiprows=2).iloc[:, :3].dropna(subset=[0]).set_index(0)
      if os.path.isfile(S1) else None)
coll = d[d['class'] == 'amplification_artefact']
print('\nDeletions reclassified as amplification artefacts:\n',
      coll[['Sequence Change', 'Details']].assign(genes=coll.genes.apply(lambda g: ', '.join(g[:8]) + (' ...' if len(g) > 8 else '')),
                                                 pops=coll.pops.apply(','.join)).to_string())
rep_rows = []
for gl, lofb in [('LoF incl. amplification artefacts', {b for bs in pd.concat([lof, coll]).bnums for b in bs}),
                 ('LoF (artefacts excluded)', set(gt.bnumber))]:
    srcs = [] if s1 is None else [('Table S1 "Glucose" column (as published; 0 for genes without data)', pd.to_numeric(s1[2], errors='coerce'))]
    for label, vals in srcs + [('repo genesfitness_Glucose.csv (original)', orig.fitval),
                        ('corrected mean IT003/IT004 (this analysis)', glc)]:
        v = vals.reindex(sorted(lofb)).dropna()
        c = pd.Series(classify(v), index=v.index).value_counts()
        rep_rows.append({'gene set': gl, 'fitness source': label, 'LoF genes with a value': len(v), **{k: int(c.get(k, 0)) for k in CATS},
                         'mean-effect_%': round(100 * c.get('mean-effect', 0) / len(v), 1),
                         'of which exactly 0 (no data)': int((v == 0).sum())})
rep_rows.append({'fitness source': 'manuscript (l. 560-566)', 'LoF genes with a value': 459, 'essential': 13, 'important': 12,
                 'mean-effect': 423, 'enhancing': 11, 'mean-effect_%': 92.2})
rep_df = pd.DataFrame(rep_rows)
rep_df.to_csv(f'{OUT}/a09{SUF}_reproduce_manuscript_counts.csv', index=False)
print('\nReproduction of the manuscript counts:\n', rep_df.to_string(index=False))

# ---------------------------------------------------------------------------------------------------
# 7. figure
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for p, g in cinfo.groupby('population'):
    ax[0].plot(g.gen, g.n_mut, 'o', ms=3, label=p)
ax[0].set_xlabel('generation'); ax[0].set_ylabel('mutations per clone'); ax[0].legend(frameon=False, fontsize=7)
ax[0].set_title('No clone shows a hypermutator jump', fontsize=9)
r0 = res[(res.set == 'all LoF') & (res.categories == 'fitness only')].iloc[0]
xs = np.arange(len(CATS))
ax[1].bar(xs - 0.2, [r0[f'{k}_background_%'] for k in CATS], 0.4, color='#bbbbbb', label='all genes (glucose)')
ax[1].bar(xs + 0.2, [r0[f'{k}_%'] for k in CATS], 0.4, color=[COLORS[k] for k in CATS], label='LoF genes in LTEE')
ax[1].set_xticks(xs, CATS, fontsize=8); ax[1].set_ylabel('% of genes'); ax[1].legend(frameon=False, fontsize=7)
ax[1].set_title(f'LoF genes (n = {int(r0.genes_in_universe)}) vs background', fontsize=9)
bins = np.linspace(-4, 2, 61)
ax[2].hist(glc.reindex(universe).clip(-4, 2), bins, density=True, color='#bbbbbb', label='all genes')
ax[2].hist(glc.reindex(u.bnumber).clip(-4, 2), bins, density=True, histtype='step', color='k', lw=1.5, label='LoF genes')
ax[2].axvline(-0.5, ls=':', c='grey'); ax[2].axvline(0.5, ls=':', c='grey')
ax[2].set_xlabel('glucose fitness'); ax[2].set_ylabel('density'); ax[2].legend(frameon=False, fontsize=7)
plt.tight_layout(); plt.savefig(f'{OUT}/a09{SUF}_ltee.png', dpi=200); plt.savefig(f'{OUT}/a09{SUF}_ltee.svg')
print('\nSaved a09_* to', OUT)
