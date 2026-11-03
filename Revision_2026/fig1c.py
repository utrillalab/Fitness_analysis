import pandas as pd, numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('FA_REPO', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # clone of utrillalab/Fitness_analysis
os.chdir(os.path.join(REPO, 'Code_Notebooks'))
# cell 12
table = pd.read_excel('../Input_Data/41587_2016_BFnbt3418_MOESM18_ESM_Mod.xlsx', sheet_name='Table S6', skiprows=[0,1])
pc = table.iloc[:, [0,2,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,74,75,76,77]]
protein_mass = table.iloc[:, [0,2,5,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,74,75,76,77]]
protein_mass.columns = pc.columns
protein_mass = protein_mass.drop(labels=2123, axis=0)
# cell 13
ft = pd.read_csv('../Input_Data/fit_organism_Keio.tsv', sep='\t')
fa = pd.DataFrame(); fa['bnumber']=ft.sysName; fa.index=fa.bnumber
for c in range(5,95,2): fa[ft.columns[c]] = ft.iloc[:,[c,c+1]].mean(axis=1).values
fa['LB'] = ft.iloc[:,[112,131]].mean(axis=1).values
sf = fa.iloc[:, [0,1,12,25,19,15,4,27,5,16,2,46]]
sf.columns = ['bnumber','Glucose','Acetate','Glucosamine','Glycerol','Pyruvate','Xylose','Mannose','Galactose','Succinate','Fructose','LB']
print('fitness genes', len(sf), '| proteomics genes', protein_mass.Bnumber.nunique())
shared = set(protein_mass.Bnumber) & set(sf.bnumber)
print('shared genes', len(shared))
def cls(f):
    if f<=-2: return 'essential'
    if -0.5<=f<=0.5: return 'no_effect'
    if -2<f<-0.5: return 'important'
    if f>0.5: return 'fitness-enhancing'
    return 'huh'
pm = protein_mass.drop_duplicates('Bnumber').set_index('Bnumber')
conds = ['LB','Glycerol','Galactose','Mannose','Glucosamine','Pyruvate','Fructose','Glucose','Xylose','Succinate','Acetate']
print('cols in proteomics?', [c for c in conds if c not in pm.columns])
rows=[]
for c in conds:
    g = sorted(shared)
    fg = pd.to_numeric(pm.loc[g, c], errors='coerce')
    cl = sf.loc[g, c].map(cls)
    tot = np.nansum(fg)
    r = {'condition': c}
    for k in ['no_effect','essential','important','fitness-enhancing','huh']:
        r[k] = 100*np.nansum(fg[cl.values==k])/tot
    r['n_genes_fg>0'] = int((fg>0).sum())
    rows.append(r)
out = pd.DataFrame(rows).set_index('condition').round(1)
print(out.to_string())
print('\nno_effect range: %.1f-%.1f' % (out.no_effect.min(), out.no_effect.max()))
print('essential range: %.1f-%.1f' % (out.essential.min(), out.essential.max()))
# huh = NaN fitness?
print('NaN fitness among shared, per cond:', {c:int(sf.loc[sorted(shared),c].isna().sum()) for c in conds})
out.to_csv(os.path.join(HERE,'fig1c_fractions.csv'))
