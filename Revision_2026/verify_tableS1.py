"""Diagnose Table S1: compare its fitness/proteomics columns with the raw data and find which replicate pairs were averaged.
Uses fig1c.py from this folder; repo clone path via FA_REPO (default: the repository root)."""
import pandas as pd, numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, 'fig1c.py')).read().split('print(\'fitness genes')[0])
pm = protein_mass.drop_duplicates('Bnumber').set_index('Bnumber')
s1 = pd.read_excel(os.path.join(HERE, 'external', 'Supp_Tables_1-19_Hidalgo_et_al.xlsx'), sheet_name='Table S1', header=None, skiprows=2)
C=['Glucose','LB','Glucosamine','Glycerol','Xylose','Succinate','Mannose','Galactose','Fructose','Acetate','Pyruvate']
s1 = s1.iloc[:, :24]; s1.columns=['bnum','gene']+['F_'+c for c in C]+['P_'+c for c in C]
s1 = s1.dropna(subset=['bnum']).set_index('bnum')
g = sorted(set(s1.index)&set(sf.index)&set(pm.index))
print('common genes', len(g))
# for each S1 fitness column, find which notebook condition it correlates best with
print('\nFITNESS: S1 column -> best-matching notebook condition (max abs diff)')
for c in C:
    a = s1.loc[g,'F_'+c].astype(float)
    best = min(['Glucose','Acetate','Glucosamine','Glycerol','Pyruvate','Xylose','Mannose','Galactose','Succinate','Fructose','LB'], key=lambda k: np.nanmean(np.abs(a.values - sf.loc[g,k].values)))
    print(f'  {c:12s} -> {best:12s} mean|diff| same-label={np.nanmean(np.abs(a.values-sf.loc[g,c].values)):.3f}')
print('\nPROTEOMICS: S1 column -> best-matching Schmidt condition')
for c in C:
    a = s1.loc[g,'P_'+c].astype(float)
    best = min(C, key=lambda k: np.nanmean(np.abs(a.values - pd.to_numeric(pm.loc[g,k],errors='coerce').values)))
    print(f'  {c:12s} -> {best}')
print('\nS1 fitness column -> closest raw column / replicate-average in fit_organism_Keio.tsv')
raw = ft.set_index('sysName')
cand = {}
for j in range(5, ft.shape[1]):
    cand[ft.columns[j]] = pd.to_numeric(raw.iloc[:, j-1] if False else ft.set_index('sysName')[ft.columns[j]], errors='coerce')
for c in ['F_'+x for x in C]:
    a = s1.loc[g, c].astype(float).values
    sc = sorted(((np.nanmean(np.abs(a - v.loc[g].values)), k) for k, v in cand.items()))[:3]
    avg = {k: v for k, v in fa.items() if k!='bnumber'}
    sa = sorted(((np.nanmean(np.abs(a - pd.to_numeric(v.loc[g]).values)), k) for k, v in avg.items()))[:2]
    print(f'  {c:12s} raw: ' + ' | '.join(f'{k} ({d:.3f})' for d,k in sc) + '  || avg: ' + ' | '.join(f'{k} ({d:.3f})' for d,k in sa))
print('\nBest pair-average match for each S1 fitness column')
names = [ft.columns[j] for j in range(5, ft.shape[1])]
M = {k: pd.to_numeric(ft.set_index('sysName')[k], errors='coerce').loc[g].values for k in names}
import itertools
for c in ['F_'+x for x in C]:
    a = s1.loc[g, c].astype(float).values
    best = min(((np.nanmax(np.abs(a - np.round((M[p]+M[q])/2, 4))), p, q) for p, q in itertools.combinations(names, 2)))
    single = min(((np.nanmax(np.abs(a - M[p])), p) for p in names))
    print(f'  {c:12s} pair maxdiff={best[0]:.4f}: {best[1]}  +  {best[2]}   | best single {single[1]} ({single[0]:.3f})')
print('\nRepo input files vs S1 and vs notebook (Acetate fitness, mean|diff|)')
fa_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'original_submission', 'Infiles')
ag = pd.read_csv(os.path.join(fa_dir,'allgen4691_fitfg_newname.csv'))
print(ag.columns.tolist()[:30])
def cls(f):
    if f!=f: return 'NA'
    if f<=-2: return 'essential'
    if -0.5<=f<=0.5: return 'no_effect'
    if -2<f<-0.5: return 'important'
    return 'fitness-enhancing'
ag = ag.set_index('bnum')
print('allgen4691 == S1 fitness?', all(np.allclose(ag.loc[g, c+'_fit'].astype(float), s1.loc[g,'F_'+c].astype(float)) for c in C))
allg = sorted(set(sf.index) & set(ag.index))
print(f'\nGenes with fitness data: {len(allg)}. Category agreement, wrong (S1/allgen4691) vs correct (notebook):')
tot_changed = 0
for c in C:
    wrong = ag.loc[allg, c+'_fit'].astype(float).map(cls); right = sf.loc[allg, c].map(cls)
    ch = (wrong != right).sum(); tot_changed += ch
    gc = pd.read_csv(os.path.join(fa_dir,'genes_categories',f'genesfitness_{c}.csv')).set_index('bnum')
    gg = sorted(set(gc.index)&set(sf.index))
    dw = np.nanmean(np.abs(gc.loc[gg,'fitval'].values - ag.loc[gg, c+'_fit'].astype(float).values))
    dr = np.nanmean(np.abs(gc.loc[gg,'fitval'].values - sf.loc[gg, c].values))
    print(f'  {c:12s} genes changing category: {ch:4d} ({100*ch/len(allg):4.1f}%) | genesfitness_{c}.csv (n={len(gc)}) matches: {"WRONG" if dw<dr else "correct"} (|d| wrong={dw:.3f}, right={dr:.3f})')
# fitness-enhancing in >=10 conditions? Table S19 is over 46 conditions -> uses fitness_averages (correct). check 108 list source
enh = pd.read_excel(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'Functional_Analysis', 'Infiles', '108_enhancing_10Conditions.xlsx'))
print('\n108_enhancing file columns:', enh.columns.tolist()[:8], len(enh))
