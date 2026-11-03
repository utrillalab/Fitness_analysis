"""Tables S17/S18 with the revised functional annotation of J. A. Farías (v3, 16 Oct 2026; R1-m12/R2-m1, Fig. 7).

Source: ../external/JoseAFarias_fig7/Supplementary-Table-functional-annotation-v3.docx (one table, 28 targets).
- "Predicted function" = the v3 "Functional Annotation" column.
- Table S17 = high confidence (HHsearch probability >= 90% and DALI Z-score >= 8), Table S18 = the rest (lower confidence and
  no hits), as stated in Methods and in the Fig. 7 legend. This moves yfaT (Z 29.3, 99.85%) to S17 and yfdY (Z 7.2, 61.8%)
  to S18; all other genes keep their table.
- Gene ID, EcoCyc description, conditions with fitness > 0.5 and t > 4, and core membership are kept from a12.
Output: ../Table_S17_S18_revised.xlsx (previous version kept as Table_S17_S18_revised_before_v3.xlsx) and
output/a32_tablesS17_S18_v3.csv
"""
import os, re, html, zipfile, shutil
import pandas as pd
from common import HERE, OUT

ROOT = os.path.join(os.path.dirname(HERE), 'external')
V3 = os.path.join(ROOT, 'JoseAFarias_fig7', 'Supplementary-Table-functional-annotation-v3.docx')
XLSX = os.path.join(os.path.dirname(HERE), 'Table_S17_S18_revised.xlsx')


def cell_text(c):
    paras = re.findall(r'<w:p[ >].*?</w:p>', c, re.S)
    return '\n'.join(''.join(html.unescape(t) for t in re.findall(r'<w:t(?: [^>]*)?>([^<]*)</w:t>', p)).strip()
                     for p in paras).strip()


x = zipfile.ZipFile(V3).read('word/document.xml').decode('utf8')
rows = []
for tbl in re.findall(r'<w:tbl>.*?</w:tbl>', x, re.S):
    for tr in re.findall(r'<w:tr[ >].*?</w:tr>', tbl, re.S):
        c = [cell_text(t) for t in re.findall(r'<w:tc>.*?</w:tc>', tr, re.S)]
        g = re.search(r'gene:?\s*(\w+)', c[1]) if len(c) > 1 else None
        if not g and len(c) > 1:  # yebG row has no "gene:"; take the name from the UniProt entry (YEBG_ECOLI -> yebG)
            u = re.search(r'\|(\w{4})_ECOLI', c[0])
            g = re.match(r'(.*)', u.group(1)[:3].lower() + u.group(1)[3]) if u else None
        if not g:
            continue
        z = re.search(r'z-score:\s*([\d.]+)', c[2], re.I)
        p = re.search(r'Probability:\s*([\d.]+)', c[3], re.I)
        rows.append({'gene': g.group(1), 'z': float(z.group(1)) if z else None, 'prob': float(p.group(1)) if p else None,
                     'annotation': c[4] if len(c) > 4 else ''})
v3 = pd.DataFrame(rows).drop_duplicates('gene').set_index('gene')
v3['table'] = ['Table S17' if (pr or 0) >= 90 and (zz or 0) >= 8 else 'Table S18' for pr, zz in zip(v3.prob, v3.z)]
assert len(v3) == 28 and (v3.table == 'Table S17').sum() == 15, v3.table.value_counts()

old = pd.concat([pd.read_excel(XLSX, sheet_name=s).assign(old_table=s) for s in ['Table S17', 'Table S18']])
old['gene'] = old.iloc[:, 1].astype(str).str.extract(r'gene\s*:?\s*(\w+)')[0].values
old = old.set_index('gene')
assert set(old.index) == set(v3.index)
new = old.join(v3[['table', 'annotation', 'prob', 'z']])
new['Predicted function'] = (new.annotation.fillna('').str.replace(r'\s*\n\s*', '; ', regex=True)
                            .replace('', 'No prediction (no hits above the retention thresholds)'))
# J. A. Farias's answers (22 Oct 2026): ydaW has two tandem winged-helix modules ("2x" in v3 = two hits);
# yibN lacks the conserved catalytic cysteine of rhodaneses
FARIAS = {'ydaW': ('Two tandem winged-helix DNA-binding modules (residues 1–54 and 78–138); 7BZH (a Sulfolobus DNA-binding '
                   'protein) aligns with each module (residues 5–53 and 82–136); the structural match to the Zα domain of '
                   'PKZ (4LB5) covers about 50 residues'),
          'yibN': 'Rhodanese-like domain (GlpE-like); putative sulfurtransferase, but no conserved catalytic cysteine'}
for g, t in FARIAS.items():
    assert g in new.index, g
    new.loc[g, 'Predicted function'] = t
moved = new[new.old_table != new.table]
print('moved:', ', '.join(f'{g} {a} -> {b}' for g, a, b in zip(moved.index, moved.old_table, moved.table)))
new.reset_index()[['gene', 'old_table', 'table', 'prob', 'z', 'Predicted function']].to_csv(
    f'{OUT}/a32_tablesS17_S18_v3.csv', index=False)

if not os.path.exists(XLSX.replace('.xlsx', '_before_v3.xlsx')):
    shutil.copy(XLSX, XLSX.replace('.xlsx', '_before_v3.xlsx'))
cols = [c for c in old.columns if c not in ('old_table',)]
with pd.ExcelWriter(XLSX) as w:
    for s in ['Table S17', 'Table S18']:
        d = new[new.table == s]
        # original order (by protein size); a moved gene keeps its position in the size-sorted v3 table
        d = d.loc[[g for g in v3.index if g in d.index]]
        d[cols].to_excel(w, sheet_name=s, index=False)
        print(s, len(d), ', '.join(d.index))
print('saved', XLSX)
