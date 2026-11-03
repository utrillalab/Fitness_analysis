"""Revised Supplementary Table S19: per-gene counts of conditions (of 46) in each fitness category, plus new
columns for the fitness-enhancing candidate set (fitness > 0.5 in >= 10 conditions, 108 genes) and the
statistically supported core (fitness > 0.5 in >= 10 of 46 conditions and t > 4 in at least one of them).
t-like statistics: Price et al. 2018 (Input_Data_Price2018/fit_t.tab); replicate t combined as
sum(t)/sqrt(n)/sqrt((n^2+n)/(2n)) (Price et al. 2018 Methods). The original workbook is not modified."""
import os
import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from common import *

ROOT = os.path.dirname(HERE)
SUPP = os.path.join(ROOT, 'external', 'Supp_Tables_1-19_Hidalgo_et_al.xlsx')
s19 = pd.read_excel(SUPP, sheet_name='Table S19')

tt = pd.read_csv(os.path.join(ROOT, 'Input_Data_Price2018', 'fit_t.tab'), sep='\t').set_index('sysName')
ft = pd.read_csv(os.path.join(REPO, 'Input_Data', 'fit_organism_Keio.tsv'), sep='\t').set_index('sysName')
fa, names46 = load_fitness_all()
t46 = {}
for col in fa.columns:
    pair = ['set2IT045 LB', 'set2IT069 LB'] if col == 'LB' else [col, ft.columns[list(ft.columns).index(col) + 1]]
    t46[col] = tt[pair].sum(axis=1) / np.sqrt(3)
t46 = pd.DataFrame(t46).reindex(fa.index)
n_enh = (fa > 0.5).sum(axis=1)
n_enh_sig = ((fa > 0.5) & (t46 > 4)).sum(axis=1)
assert (s19.set_index('bnumber').sum_of_enhancing.reindex(fa.index) == n_enh).all(), 'S19 counts differ'

s19['sum_of_enhancing_significant (t > 4)'] = s19.bnumber.map(n_enh_sig).astype(int)
s19['recurrent_enhancing_candidate (>=10)'] = np.where(s19.sum_of_enhancing >= 10, 'yes', '')
s19['significant_core (>=10, t > 4 in >=1)'] = np.where((s19.sum_of_enhancing >= 10) & (s19['sum_of_enhancing_significant (t > 4)'] >= 1), 'yes', '')
s19['possible_TnSeq_artefact'] = np.where(  # essential in Keio, PEC and DEG (R1-m13, R2-m2)
    s19.gene_name.isin(['rpsG', 'rpsR']) & (s19.sum_of_enhancing >= 10), 'yes', '')
n_cand = (s19['recurrent_enhancing_candidate (>=10)'] == 'yes').sum()
n_core = (s19['significant_core (>=10, t > 4 in >=1)'] == 'yes').sum()
print(f'candidates: {n_cand}, significant core: {n_core}, '
      f'core not in candidates: {((s19["significant_core (>=10, t > 4 in >=1)"] == "yes") & (s19["recurrent_enhancing_candidate (>=10)"] != "yes")).sum()}')

wb = Workbook(); ws = wb.active; ws.title = 'Table S19'
ws['A1'] = 'Supplementary Table S19. Number of conditions (of 46) in which each gene falls in each fitness category.'
ws['A1'].font = Font(bold=True)
ws.append(list(s19.columns))
for cell in ws[2]:
    cell.font = Font(bold=True); cell.alignment = Alignment(wrap_text=True, vertical='top')
for row in s19.itertuples(index=False):
    ws.append(list(row))
for n in ['', 'New columns (revision):',
          'sum_of_enhancing_significant: conditions with fitness > 0.5 and t > 4 (significance criterion of Price et al. 2018).',
          f'recurrent_enhancing_candidate: fitness > 0.5 in >= 10 conditions ({n_cand} genes).',
          f'significant_core: fitness > 0.5 in >= 10 of 46 conditions and t > 4 in at least one of them ({n_core} genes).',
          'possible_TnSeq_artefact: rpsG and rpsR, the only candidates essential in the Keio collection, PEC and DEG '
          '(they may appear fitness-enhancing because of insertions tolerated near the 3\' end or polar effects).']:
    ws.append([n])
for col, w in zip('ABCDEFGHIJK', [10, 12, 50, 12, 12, 12, 12, 16, 16, 16, 18]):
    ws.column_dimensions[col].width = w
ws.freeze_panes = 'D3'
ws.auto_filter.ref = f'A2:K{2 + len(s19)}'
out = os.path.join(ROOT, 'Table_S19_revised.xlsx')
wb.save(out)
print('saved', out)
print(s19[s19['significant_core (>=10, t > 4 in >=1)'] == 'yes'].sort_values('sum_of_enhancing_significant (t > 4)', ascending=False)
      [['gene_name', 'sum_of_enhancing', 'sum_of_enhancing_significant (t > 4)', 'possible_TnSeq_artefact']].head(8).to_string(index=False))
