"""Revised Supplementary Tables S2-S12 (GO enrichment per condition) and a new table for the R1-M1 sensitivity analysis.

S2-S12: ViSEAGO results of a15_viseago.R, setting 'corrected' (correct fitness values, glucosamine HCl, the 287 essential
genes without fitness data added to the essential category; GO annotation org.EcK12.eg.db 3.22.0 / GO.db 3.23.1),
in the layout of the original sheets (title row + ViSEAGO columns, "No_effect" renamed "Mean_effect").
New table: essential-category GO terms with and without the 287 genes, per condition, and the 287-gene list
(a14_essential_nofit_list.py).
Output: ../Tables_S2_S12_revised.xlsx, ../Table_SX_essential_without_fitness.xlsx
"""
import os
import pandas as pd
from common import *

ORDER = [('S2', 'Glucose'), ('S3', 'LB'), ('S4', 'Glucosamine'), ('S5', 'Glycerol'), ('S6', 'Xylose'), ('S7', 'Mannose'),
         ('S8', 'Galactose'), ('S9', 'Fructose'), ('S10', 'Acetate'), ('S11', 'Pyruvate'), ('S12', 'Succinate')]
NOTE = ('Revised: fitness = mean of the two replicates of each condition (glucosamine = D-glucosamine hydrochloride); the 287 '
        'genes without fitness data annotated as essential (PEC, DEG or Keio) are added to the essential category; '
        'GO annotation org.EcK12.eg.db 3.22.0, GO.db 3.23.1 (Bioconductor 3.23); topGO elim + Fisher, p < 0.01.')
base = os.path.dirname(HERE)
with pd.ExcelWriter(os.path.join(base, 'Tables_S2_S12_revised.xlsx')) as w:
    for tab, cond in ORDER:
        d = pd.read_csv(os.path.join(OUT, 'viseago', 'corrected', f'BP_Results_{cond}.xls'), sep='\t')
        d.columns = [c.replace('No_effect', 'Mean_effect') for c in d.columns]
        title = pd.DataFrame([[f'Supplementary Table {tab}. GO terms enriched in the {cond} condition. ' + NOTE]])
        title.to_excel(w, sheet_name=f'Table {tab}', index=False, header=False)
        d.to_excel(w, sheet_name=f'Table {tab}', index=False, startrow=1)
        print(f'Table {tab} ({cond}): {len(d)} GO terms')

e = pd.read_csv(os.path.join(OUT, 'a15b_essential_terms_with_vs_without_287.csv'))
e['enriched_with_287'] = e.p_with_287.notna(); e['enriched_without_287'] = e.p_without_287.notna()
summ = (e.groupby(['GO.ID', 'term']).agg(conditions_with_287=('enriched_with_287', 'sum'),
                                         conditions_without_287=('enriched_without_287', 'sum'))
         .reset_index().sort_values(['conditions_with_287', 'conditions_without_287'], ascending=False))
genes = pd.read_csv(os.path.join(OUT, 'a14_essential_from_nofit_list.csv'))
with pd.ExcelWriter(os.path.join(base, 'Table_SX_essential_without_fitness.xlsx')) as w:
    pd.DataFrame([['Genes without RB-TnSeq fitness data that are annotated as essential in at least one of three sources '
                   '(PEC, DEG, Keio collection) and were added to the essential category for GO enrichment.']]
                 ).to_excel(w, sheet_name='Genes (287)', index=False, header=False)
    genes.to_excel(w, sheet_name='Genes (287)', index=False, startrow=1)
    pd.DataFrame([['Essential-category GO terms (p < 0.01) in the 11 conditions with and without the 287 genes '
                   '(sensitivity analysis, Reviewer 1 major comment 1).']]).to_excel(w, sheet_name='GO with vs without', index=False, header=False)
    summ.to_excel(w, sheet_name='GO with vs without', index=False, startrow=1)
print(f'Table SX: {len(genes)} genes; {len(summ)} essential GO terms '
      f'({(summ.conditions_without_287 == 0).sum()} only with the 287 genes, {(summ.conditions_with_287 == 0).sum()} only without)')
