"""Write the corrected Supplementary Table S1 in the layout of the original sheet (R2-M1, R1-M1).

Input: regenerated_inputs/Table_S1_corrected.csv (from regenerate_inputs.py).
Output: Table_S1_revised.xlsx — sheet "Table S1" (same two header rows and column order as the original; genes without
RB-TnSeq data are written as "NA" instead of 0) and sheet "Notes" (what changed and how Fig. 1C is computed).
Also checks that the Fig. 1C proteome fractions can be recomputed from this table alone.
"""
import os
import numpy as np, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment

HERE = os.path.dirname(os.path.abspath(__file__))
t = pd.read_csv(os.path.join(HERE, 'regenerated_inputs', 'Table_S1_corrected.csv'))
C = ['Glucose', 'LB', 'Glucosamine', 'Glycerol', 'Xylose', 'Succinate', 'Mannose', 'Galactose', 'Fructose', 'Acetate', 'Pyruvate']

# ---- check: Fig. 1C from this table alone (genes with fitness and proteomics; mass-weighted fractions) ----------
def cls(f):
    return np.select([f <= -2, f < -0.5, f <= 0.5, f > 0.5], ['Essential', 'Important', 'Mean-effect', 'Enhancing'], 'NA')

rows = []
for c in ['LB', 'Glycerol', 'Galactose', 'Mannose', 'Glucosamine', 'Pyruvate', 'Fructose', 'Glucose', 'Xylose', 'Succinate', 'Acetate']:
    f = t[f'Fitness_{c}']; m = pd.to_numeric(t[f'Proteome_fg_{c}'], errors='coerce')
    ok = f.notna() & m.notna(); k = cls(f[ok].values); tot = m[ok].sum()
    rows.append({'Condition': c, **{x: round(100 * m[ok][k == x].sum() / tot, 1) for x in ['Mean-effect', 'Essential', 'Important', 'Enhancing']}})
fig1c = pd.DataFrame(rows)
print(fig1c.to_string(index=False))

# ---- workbook ----------------------------------------------------------------------------------------------------
wb = Workbook()
ws = wb.active; ws.title = 'Table S1'
ws.append(['bnum', 'gene'] + C + C)
ws.append([None, None, 'Fitness values (mean of the two replicates of each condition; NA = no RB-TnSeq data)'] + [None] * 10 +
          ['Proteomics (protein mass, fg per cell; Schmidt et al. 2016)'] + [None] * 10)
for r in t.itertuples(index=False):
    d = r._asdict()
    fit = ['NA' if pd.isna(d[f'Fitness_{c}']) else round(float(d[f'Fitness_{c}']), 4) for c in C]
    prot = [None if pd.isna(d[f'Proteome_fg_{c}']) else float(d[f'Proteome_fg_{c}']) for c in C]
    ws.append([d['bnum'], d['gene']] + fit + prot)
for cell in ws[1] + ws[2]:
    cell.font = Font(bold=True)
ws.freeze_panes = 'C3'

notes = wb.create_sheet('Notes')
lines = [
    ['Supplementary Table S1 (revised). Integration of fitness and proteomics data.'],
    [],
    ['Fitness', 'RB-TnSeq gene fitness from Price et al. (2018), E. coli BW25113 (Keio library). Each value is the mean of the two replicate experiments of the same condition.'],
    ['', 'Glucosamine = D-glucosamine hydrochloride (experiments IT063/IT064), the condition used for proteomics by Schmidt et al. (2016).'],
    ['', f'NA = gene without RB-TnSeq fitness data ({int(t.filter(like="Fitness_").isna().all(axis=1).sum())} genes). In the previous version these genes were coded as 0.'],
    ['Proteomics', 'Protein mass per cell (fg) from Schmidt et al. (2016), Table S6; unchanged from the previous version.'],
    ['Categories', 'Essential: fitness <= -2; Important: -2 < fitness < -0.5; Mean-effect: -0.5 <= fitness <= 0.5; Fitness-enhancing: fitness > 0.5.'],
    ['Fig. 1C', 'Proteome fractions are computed over genes with both fitness and proteomics data in each condition (genes with NA fitness are excluded), as the protein mass of each category divided by the total mass of those genes.'],
    [],
    ['Fig. 1C proteome fractions recomputed from this table (%)'],
    list(fig1c.columns),
] + fig1c.values.tolist()
for l in lines:
    notes.append(l)
notes['A1'].font = Font(bold=True); notes.cell(row=10, column=1).font = Font(bold=True)
for cell in notes[11]:
    cell.font = Font(bold=True)
notes.column_dimensions['A'].width = 14; notes.column_dimensions['B'].width = 120
for row in notes.iter_rows(min_row=3, max_row=8):
    row[1].alignment = Alignment(wrap_text=True, vertical='top')

out = os.path.join(HERE, 'Table_S1_revised.xlsx')
wb.save(out)
print('saved', out, '| genes:', len(t))
