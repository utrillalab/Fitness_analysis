"""Reproducible version of the 287-gene list used to extend the essential category in the GO enrichment
(Lizeth's notebook Fig_suppl01_intersection_essential_genes.py, which hard-codes server paths and
cannot be run top to bottom).

Rule (as in the notebook): a gene without RB-TnSeq fitness data is called essential if it is annotated
as essential in at least one of PEC, DEG or the 276 Keio essential genes (gene names as in RegulonDB 12).
The no-fitness set is rebuilt from the revised Table S1 (genes with NA fitness). The output adds one
column per source so that the Methods can report how many genes each source contributes."""
import os
import pandas as pd
from common import *

ROOT = os.path.dirname(HERE)
TOP = os.path.join(ROOT, 'external', 'essential_lists')
pec = set(pd.read_csv(os.path.join(TOP, 'genes_pec_synon.csv')).GenesReg12.dropna().astype(str))
deg = set(pd.read_csv(os.path.join(TOP, 'DEG_synonym.csv')).Genes_R12.dropna().astype(str))
keio = set(pd.read_csv(os.path.join(TOP, 'genes_keio276_synon.csv')).GenesReg12.dropna().astype(str))
gl = pd.read_excel(os.path.join(TOP, '41589_2020_593_MOESM2_ESM.xlsx'), sheet_name='Table 1', header=0)
gl_core = set(gl.loc[gl.Classification == 'Core Essential', 'Name'].astype(str))
print(f'PEC {len(pec)}, DEG {len(deg)}, Keio essential {len(keio)}, Lastiri-Pancardo 2020 core {len(gl_core)}')

# genes without fitness data = genes of Table S1 whose fitness is NA in all conditions
s1 = pd.read_excel(os.path.join(ROOT, 'Table_S1_revised.xlsx'), sheet_name=0, header=None)
hdr = next(i for i in range(10) if str(s1.iat[i, 0]).strip().lower() in ('bnum', 'bnumber'))
s1 = s1.iloc[hdr + 1:].reset_index(drop=True)
s1 = s1[s1[0].astype(str).str.match(r'^b\d{4}$')]
fitcols = list(range(2, 13))
nofit = s1[s1[fitcols].apply(pd.to_numeric, errors='coerce').isna().all(axis=1)]
print('genes without fitness data in revised Table S1:', len(nofit))

# names: Table S1 gene name plus EcoCyc name, to match RegulonDB-12 names
eco = pd.read_csv(os.path.join(REPO, 'Input_Data', 'All-genes-of-E.-coli-K-12-substr-MG1655.txt'), sep='\t')
eco_name = eco.dropna(subset=['Accession-1']).drop_duplicates('Accession-1').set_index('Accession-1')['Gene Name']
def names_for(b, g):
    n = {str(g), str(eco_name.get(b, g))}
    return {x for x in n if x and x != 'nan'}
rows = []
for b, g in zip(nofit[0], nofit[1]):
    n = names_for(b, g)
    rows.append({'bnumber': b, 'gene': g, 'PEC': bool(n & pec), 'DEG': bool(n & deg), 'Keio_essential': bool(n & keio),
                 'Lastiri2020_core (glucose only, not used)': bool(n & gl_core)})
t = pd.DataFrame(rows)
t['essential_in_>=1_source'] = t[['PEC', 'DEG', 'Keio_essential']].any(axis=1)
sel = t[t['essential_in_>=1_source']]
print('reproduced list:', len(sel))

orig = set(pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'essential_287gdeNofit_otherConds.csv')).gene.astype(str))
mine = set(sel.gene.astype(str)) | {str(eco_name.get(b)) for b in sel.bnumber}
print('original list genes not reproduced:', sorted(orig - mine))
print('reproduced genes not in original list:', sorted(set(sel.gene.astype(str)) - orig))

print('\nContribution of each source (genes among the selected):')
print(sel[['PEC', 'DEG', 'Keio_essential']].sum().to_string())
combo = sel[['PEC', 'DEG', 'Keio_essential']].apply(lambda r: '+'.join(k for k, v in r.items() if v), axis=1)
print(combo.value_counts().to_string())

k = pd.read_excel(os.path.join(REPO, 'Input_Data', 'Map_KeioCollection.xlsx'), sheet_name='Data2')
kb = set(k['Escherichia coli MG1655 B id'].dropna().astype(str))
sel = sel.assign(has_Keio_mutant=sel.bnumber.isin(kb), sources=combo)
print('\nSelected genes with a Keio mutant, and their sources:')
print(sel[sel.has_Keio_mutant][['gene', 'sources']].to_string(index=False))

t.to_csv(f'{OUT}/a14_nofit_genes_essential_sources.csv', index=False)
sel.to_csv(f'{OUT}/a14_essential_from_nofit_list.csv', index=False)
