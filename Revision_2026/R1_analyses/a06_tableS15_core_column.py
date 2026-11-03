"""Revised Supplementary Table S15: GO terms enriched among the 108 candidate fitness-enhancing genes,
with new columns marking the genes of each term that belong to the statistically supported core
(genes with fitness > 0.5 in >= 10 of 46 conditions and t > 4 in at least one of them; see a05_t_statistic.py).
Rows are sorted by -log10 p-value (R1-m11). GO terms: enrichment rerun with the current GO annotation (a30_table_S15_GO.R,
same parameters as the original analysis), as for Tables S2-S12 and S14. The original supplementary workbook is not modified."""
import os
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter
from common import *

s15 = pd.read_csv(f'{OUT}/viseago/enhancing_108/BP_Results_enhancing_108.xls', sep='	')
cand = pd.read_csv(f'{OUT}/a05_enhancing_108_candidates_32_core.csv')

# symbols for each b-number: fitness-table name + EcoCyc name (GO tables use current symbols, e.g. rsmG for gidB)
eco = pd.read_csv(os.path.join(REPO, 'Input_Data', 'All-genes-of-E.-coli-K-12-substr-MG1655.txt'), sep='\t')
eco_name = eco.dropna(subset=['Accession-1']).drop_duplicates('Accession-1').set_index('Accession-1')['Gene Name']
sym2b = {}
for b, g in zip(cand.bnumber, cand.gene):
    for s in {str(g), str(eco_name.get(b, g))}:
        sym2b[s.lower()] = b
core_b = set(cand.loc[cand.significant_core, 'bnumber'])
artefact_b = set(cand.loc[cand.gene.isin(['rpsG', 'rpsR']), 'bnumber'])  # essential in Keio, PEC and DEG (R1-m13, R2-m2)

core_cnt, core_sym, cand_only, flagged, unmatched = [], [], [], [], set()
for genes in s15['Fitness_enhancing.Significant_genes_symbol'].astype(str):
    syms = [s.strip() for s in genes.split(';') if s.strip()]
    b = [sym2b.get(s.lower()) for s in syms]
    unmatched |= {s for s, x in zip(syms, b) if x is None}
    in_core = [s for s, x in zip(syms, b) if x in core_b]
    core_cnt.append(f'{len(in_core)}/{len(syms)}')
    core_sym.append(';'.join(in_core))
    cand_only.append(';'.join(s for s, x in zip(syms, b) if x not in core_b))
    flagged.append(';'.join(s for s, x in zip(syms, b) if x in artefact_b))
s15['Core_genes (core/term genes)'] = core_cnt
s15['Core_genes_symbol'] = core_sym
s15['Candidate_only_genes_symbol'] = cand_only
s15['Possible_TnSeq_artefact'] = flagged
s15 = s15.sort_values('Fitness_enhancing.-log10_pvalue', ascending=False).reset_index(drop=True)
if unmatched:
    print('WARNING - symbols not matched to the 108 candidates:', sorted(unmatched))

# write workbook: title row, header, data, notes
wb = Workbook(); ws = wb.active; ws.title = 'Table S15'
ws['A1'] = ('Supplementary Table S15. GO terms enriched among the 108 candidate fitness-enhancing genes '
            '(fitness > 0.5 in >= 10 of 46 conditions; topGO elim + Fisher, p < 0.01, current GO annotation), sorted by -log10 p-value.')
ws['A1'].font = Font(bold=True)
ws.append(list(s15.columns))
for cell in ws[2]:
    cell.font = Font(bold=True); cell.alignment = Alignment(wrap_text=True, vertical='top')
for row in s15.itertuples(index=False):
    ws.append(list(row))
notes = [
    '',
    'New columns (revision):',
    'Core_genes: genes of the term that belong to the statistically supported core, i.e. fitness > 0.5 in >= 10 of 46 conditions and t > 4 in at least one of them, '
    'using the per-gene t-like statistics and significance criterion of Price et al. (2018).',
    'Candidate_only_genes: genes of the term that pass fitness > 0.5 in >= 10 conditions but never reach t > 4.',
    'Possible_TnSeq_artefact: rpsG and rpsR, the only candidates essential in the Keio collection, PEC and DEG; they may '
    'appear fitness-enhancing because of insertions tolerated near the 3\' end or polar effects (see response to Reviewer 2).',
]
for n in notes:
    ws.append([n])
widths = {'A': 12, 'B': 38, 'C': 60, 'D': 16, 'E': 12, 'F': 12, 'G': 30, 'H': 45, 'I': 14, 'J': 30, 'K': 45, 'L': 18}
for col, w in widths.items():
    ws.column_dimensions[col].width = w
for r in ws.iter_rows(min_row=3, max_row=2 + len(s15)):
    for cell in r:
        cell.alignment = Alignment(wrap_text=True, vertical='top')
ws.freeze_panes = 'C3'
out = os.path.join(os.path.dirname(HERE), 'Table_S15_revised.xlsx')
wb.save(out)
s15.to_csv(f'{OUT}/a06_table_S15_revised.csv', index=False)
print(s15[['GO.ID', 'term', 'Fitness_enhancing.-log10_pvalue', 'Core_genes (core/term genes)', 'Core_genes_symbol',
           'Possible_TnSeq_artefact']].to_string())
print('saved', out)
