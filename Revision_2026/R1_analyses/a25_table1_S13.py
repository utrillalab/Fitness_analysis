"""Revised Table 1 and Supplementary Table S13 (glucose), from the corrected GO enrichment (a15, setting 'corrected').

Method as in Functional_Analysis/Fig04_ViseagoResult_glucose_top15.ipynb, with the corrections of the revision:
  - GO terms: output/viseago/corrected/BP_Results_Glucose.xls (p < 0.01 in a category);
  - mass of a term = summed protein mass of its genes in that category / total measured proteome (Schmidt et al. 2016,
    glucose), without dividing by the number of GO annotations (R1-m7); the same function as Fig. 4 (a20_fig4.terms_for);
  - Table 1 column "proteome mass" counts each gene once per category (union of the genes of its enriched terms).
Table S13 keeps the original layout: the 15 terms with the largest mass per category (all terms if fewer), with a
"Functional Category". Categories are taken from the original notebook mapping (cell 17); terms not in that mapping get a
category proposed from the term name (column "Category source" = "proposed") for review by L. Soto-Avila.
A second sheet lists every enriched glucose term with its mass.
Output: ../Table_1_revised.xlsx, ../Table_S13_revised.xlsx
"""
import os, re
import numpy as np, pandas as pd
from common import *
import a20_fig4 as f4

LIZ = {
    "Metabolic process": ['GO:0006355', 'GO:0006090', 'GO:0009185', 'GO:0006760', 'GO:0006310', 'GO:0022900', 'GO:0017148',
                          'GO:0006575', 'GO:0006163', 'GO:0046496'],
    "Cellular process": ['GO:0016226', 'GO:0051301', 'GO:0032506', 'GO:0042710', 'GO:0019058'],
    "Catabolic process": ['GO:0044282', 'GO:1901565', 'GO:0046395', 'GO:0044275', 'GO:0006402'],
    "Biosynthetic process": ['GO:1901566', 'GO:0043604', 'GO:0008652', 'GO:1901607', 'GO:0044272', 'GO:0000097',
                             'GO:0009067', 'GO:0009070', 'GO:0019438', 'GO:0046394', 'GO:0009244', 'GO:0072528',
                             'GO:1901576', 'GO:0044249', 'GO:0044283', 'GO:0042364'],
    "Translation": ['GO:0006412', 'GO:0000028', 'GO:0000027'],
    "Transport": ['GO:0055085', 'GO:0015849', 'GO:0006865', 'GO:0046942', 'GO:1902600', 'GO:0008643', 'GO:0009401',
                  'GO:0034219', 'GO:0030001', 'GO:0098660', 'GO:0015833', 'GO:0098655', 'GO:0098739', 'GO:1905039',
                  'GO:0000041'],
    "Response to stimulus": ['GO:0009628', 'GO:0009314', 'GO:0006974', 'GO:1901700', 'GO:0062197'],
    "Biological regulation": ['GO:0050790', 'GO:0050789'],
    "Interaction": ['GO:0044403', 'GO:0044010', 'GO:0051701'],
}
GO2CAT = {g: c for c, gs in LIZ.items() for g in gs}
RULES = [(r'translation|ribosom|protein-RNA complex|tRNA aminoacylation', 'Translation'), (r'detoxification', 'Response to stimulus'), (r'transport|import|export|secretion', 'Transport'),
         (r'biosynthetic', 'Biosynthetic process'), (r'catabolic', 'Catabolic process'),
         (r'response to|DNA repair|DNA damage', 'Response to stimulus'), (r'regulation', 'Biological regulation'),
         (r'interaction with|symbio', 'Interaction'), (r'metabolic process', 'Metabolic process')]


def propose(term):
    for pat, cat in RULES:
        if re.search(pat, term, re.I):
            return cat
    return 'Cellular process'


LABEL = {'essential': 'Essential', 'important': 'Important', 'mean-effect': 'Mean-effect', 'enhancing': 'Fitness-enhancing'}
t = f4.terms_for('Glucose')
mass = f4.mass['Glucose'].dropna(); tot = mass.sum()

# genes per category used for the GO enrichment (same input as a15 'corrected')
inp = pd.read_csv(os.path.join(os.path.dirname(HERE), 'regenerated_inputs', 'genes_categories', 'genesfitness_Glucose.csv'))
n287 = pd.read_csv(os.path.join(REPO, 'Functional_Analysis', 'Infiles', 'essential_287gdeNofit_otherConds.csv'))
cat_genes = {'essential': len(inp[inp.Categ == 'Essen']) + len(n287), 'important': (inp.Categ == 'Important').sum(),
             'mean-effect': inp.Categ.isin(['No_effect', 'No_effect2']).sum(), 'enhancing': (inp.Categ == 'Enhancing').sum()}

rows, allg = [], set()
for k in CATS:
    s = t[t.category == k]
    g = set().union(*s.genes) if len(s) else set()
    allg |= g
    rows.append({'Fitness category': LABEL[k], 'Total genes': cat_genes[k], 'Enriched GO terms': len(s),
                 'Genes within enriched GO terms': len(g), 'Proteins with data': int(mass.index.isin(list(g)).sum()),
                 'Proteome mass (%) across all enriched GO terms': round(100 * mass.reindex(list(g)).dropna().sum() / tot, 2)})
tab1 = pd.DataFrame(rows)
n_terms = t['GO.ID'].nunique(); shared = len(t) - n_terms
tab1.loc[len(tab1)] = {'Fitness category': 'Total', 'Total genes': sum(cat_genes.values()), 'Enriched GO terms': n_terms,
                       'Genes within enriched GO terms': len(allg), 'Proteins with data': int(mass.index.isin(list(allg)).sum()),
                       'Proteome mass (%) across all enriched GO terms': round(tab1.iloc[:, -1].sum(), 2)}
shared_terms = t[t.duplicated('GO.ID', keep=False)].groupby(['GO.ID', 'term']).category.apply(', '.join).reset_index()

note = (f'Table 1 (revised). Glucose. Total genes: genes of each category used for the GO enrichment (essential includes the '
        f'287 genes without fitness data annotated as essential; the remaining {4690 - sum(cat_genes.values())} genes without '
        f'fitness data are not assigned). Enriched GO terms: topGO elim + Fisher, p < 0.01; {shared} term(s) enriched in two '
        f'categories, so the total counts unique terms. Genes within enriched GO terms: union of the genes of the enriched '
        f'terms of the category. Proteins with data: of those genes, with proteomics data (Schmidt et al. 2016). Proteome '
        f'mass: summed mass of those proteins, each counted once, as % of the measured proteome.')
base = os.path.dirname(HERE)
with pd.ExcelWriter(os.path.join(base, 'Table_1_revised.xlsx')) as w:
    pd.DataFrame([[note]]).to_excel(w, sheet_name='Table 1', index=False, header=False)
    tab1.to_excel(w, sheet_name='Table 1', index=False, startrow=2)
    shared_terms.to_excel(w, sheet_name='Shared terms', index=False)

# ---- Table S13 -----------------------------------------------------------------------------------------
t = t.assign(**{'Fitness category': t.category.map(LABEL), 'Proteome Fraction (Glucose)': t['mass_%'].round(3)})
t['Functional Category'] = [GO2CAT.get(g) or propose(term) for g, term in zip(t['GO.ID'], t.term)]
t['Category source'] = np.where(t['GO.ID'].isin(GO2CAT), 'original', 'proposed')
cols = ['Fitness category', 'GO.ID', 'Proteome Fraction (Glucose)', 'term', 'Functional Category',
        'n_genes', 'n_with_proteomics']  # 'Category source' only for the printout (categories reviewed by L. Soto-Avila)
top = pd.concat([t[t.category == k].sort_values('mass_%', ascending=False).head(15) for k in CATS])[cols]
allt = pd.concat([t[t.category == k].sort_values('mass_%', ascending=False) for k in CATS])[cols[:4] + cols[5:]]
with pd.ExcelWriter(os.path.join(base, 'Table_S13_revised.xlsx')) as w:
    pd.DataFrame([['Supplementary Table S13 (revised). Enriched GO processes in glucose with the largest proteome mass '
                   '(top 15 per fitness category), each assigned to a single primary biological process. Proteome '
                   'fraction: summed mass of the genes of the term in that category, % of the measured proteome.']]
                 ).to_excel(w, sheet_name='Table S13', index=False, header=False)
    top.rename(columns={'GO.ID': 'GO', 'n_genes': 'Genes', 'n_with_proteomics': 'Genes with proteomics'}).to_excel(
        w, sheet_name='Table S13', index=False, startrow=2)
    allt.rename(columns={'GO.ID': 'GO'}).to_excel(w, sheet_name='All glucose terms', index=False)
print(tab1.to_string(index=False)); print('shared terms:\n', shared_terms.to_string(index=False))
print('\nS13:', len(top), 'rows;', (t.loc[top.index, 'Category source'] == 'proposed').sum(), 'categories assigned in the revision')
