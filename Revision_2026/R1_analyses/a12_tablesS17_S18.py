"""Tables S17/S18 (structure-based function predictions for uncharacterized fitness-enhancing genes), revised.

- Fixes the Gene ID of ycbZ (listed as b1758, which is ynjF; ycbZ = b0955).
- Adds, as in Tables S15/S19, the number of conditions (of 46) with fitness > 0.5 and t > 4, and whether the gene is in
  the statistically supported core (fitness > 0.5 in >= 10 conditions and t > 4 in >= 1; a05_t_statistic.py). All genes are candidates (108).
Output: ../Table_S17_S18_revised.xlsx and output/a12_tablesS17_S18_check.csv
"""
import os, re
import pandas as pd
from common import *

SUPP = os.path.join(os.path.dirname(HERE), 'external', 'Supp_Tables_1-19_Hidalgo_et_al.xlsx')
eco = pd.read_csv(os.path.join(REPO, 'Input_Data', 'All-genes-of-E.-coli-K-12-substr-MG1655.txt'), sep='\t')
n2b = dict(zip(eco['Gene Name'].astype(str), eco['Accession-1']))
core = pd.read_csv(os.path.join(OUT, 'a05_enhancing_108_candidates_32_core.csv')).set_index('bnumber')

CORE = f'Statistically supported core ({int(core.significant_core.sum())} genes)'
checks, out = [], {}
for sheet, new_name in [('Table_S17', 'Table S17'), ('Table_S18', 'Table S18')]:
    d = pd.read_excel(SUPP, sheet_name=sheet)
    d.columns = [c.strip() for c in d.columns]
    gene = d.iloc[:, 1].astype(str).str.extract(r'gene\s*:?\s*([A-Za-z0-9_]+)')[0]
    correct = gene.map(n2b)
    for g, given, b in zip(gene, d['Gene ID'].astype(str).str.strip(), correct):
        checks.append({'table': new_name, 'gene': g, 'Gene ID given': given, 'Gene ID (EcoCyc)': b, 'changed': given != b})
    d['Gene ID'] = correct.values
    d['Conditions with fitness > 0.5 and t > 4 (of 46)'] = [int(core.at[b, 'n_conditions_fitness>0.5_and_t>4']) for b in correct]
    d[CORE] = ['yes' if core.at[b, 'significant_core'] else 'no' for b in correct]
    out[new_name] = d
chk = pd.DataFrame(checks)
chk.to_csv(os.path.join(OUT, 'a12_tablesS17_S18_check.csv'), index=False)
print(chk[chk.changed].to_string(index=False))
for k, d in out.items():
    print(k, len(d), 'genes;', (d[CORE] == 'yes').sum(), 'in the core:',
          ', '.join(d.loc[d[CORE] == 'yes'].iloc[:, 1].str.extract(r'gene\s*:?\s*(\w+)')[0]))
with pd.ExcelWriter(os.path.join(os.path.dirname(HERE), 'Table_S17_S18_revised.xlsx')) as w:
    for k, d in out.items():
        d.to_excel(w, sheet_name=k, index=False)
print('saved Table_S17_S18_revised.xlsx')
