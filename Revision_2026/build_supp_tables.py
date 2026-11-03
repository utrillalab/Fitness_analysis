"""Revised Supplementary Tables (one workbook), numbering decided on 14 Oct 2026 ("minimum change"):
S1-S20 keep their numbers (revised versions where they exist; S20 = statistical tests of the revision), and the new
tables are appended in order of first citation:
  S21  287 genes without fitness data annotated as essential (R1-M1)        Table_SX_essential_without_fitness.xlsx
  S22  GO terms of the nine-species analysis grouped into themes (R1-m16)  R1_analyses/output/a28_fig_S6_terms_themes.csv
  S23  168 highly expressed mean-effect genes of P. putida (R1-M7)          Pp_FigureS5_share/results/FigureS5_tableC.tsv
  S24  COG enrichment of those genes (R1-M7)                               Pp_FigureS5_share/results/TableSX_COG_enrichment.tsv
  S25  genes with loss-of-function mutations in the 12 LTEE populations     R1_analyses/output/a09_all_ltee_lof_genes.csv
       (R1-M6, R2-M3; enhancing candidate/core columns updated to the 108 candidates / 71 core genes of option B)
Values are copied (no formulas); a title row and a ReadMe sheet describe each table.
Output: Revision_2026/Supplementary_Tables_revised.xlsx
"""
import os
import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
R = lambda *p: os.path.join(HERE, *p)
OUT = R('Supplementary_Tables_revised.xlsx')
orig = openpyxl.load_workbook(os.path.join(HERE, 'external', 'Supp_Tables_1-19_Hidalgo_et_al.xlsx'), read_only=True)

wb = openpyxl.Workbook()
readme = wb.active; readme.title = 'ReadMe'
readme.append(['Table', 'Description']); readme['A1'].font = readme['B1'].font = Font(bold=True)
BOLD = Font(bold=True)


def rows_of(path, sheet):
    w = openpyxl.load_workbook(path, read_only=True)
    return [list(r) for r in w[sheet].iter_rows(values_only=True)]


def add(name, rows, desc, title=None, header_row=None):
    """New sheet with rows (values); title is an optional first line; header_row (1-based, after the title) is bold."""
    ws = wb.create_sheet(name[:31])
    if title:
        ws.append([title]); ws['A1'].font = BOLD
    for r in rows:
        ws.append([None if (isinstance(v, float) and pd.isna(v)) else v for v in r])
    h = (header_row or 1) + (1 if title else 0)
    for c in ws[h]:
        c.font = BOLD; c.alignment = Alignment(wrap_text=True, vertical='top')
    ws.freeze_panes = ws.cell(row=h + 1, column=1)
    for col in ws.columns:
        L = max((len(str(c.value)) for c in list(col)[h - 1:h + 50] if c.value is not None), default=8)
        ws.column_dimensions[col[0].column_letter].width = min(max(10, L + 2), 60)
    readme.append([name.split(' (')[0], desc])
    return ws


def df_rows(df):
    return [list(df.columns)] + df.astype(object).where(df.notna(), None).values.tolist()


# S1-S19 ---------------------------------------------------------------------------------------------------------
add('Table S1', rows_of(R('Table_S1_revised.xlsx'), 'Table S1'),
    'Integration of fitness and proteomics data (revised: genes without fitness data coded as NA; correct replicate '
    'pairing; glucosamine = D-glucosamine hydrochloride).')
add('Table S1 notes', rows_of(R('Table_S1_revised.xlsx'), 'Notes'), 'Notes to Table S1 (columns and Fig. 1C calculation).')
cond = {2: 'Glucose', 3: 'LB', 4: 'Glucosamine', 5: 'Glycerol', 6: 'Xylose', 7: 'Mannose', 8: 'Galactose',
        9: 'Fructose', 10: 'Acetate', 11: 'Pyruvate', 12: 'Succinate'}
for k, c in cond.items():
    add(f'Table S{k}', rows_of(R('Tables_S2_S12_revised.xlsx'), f'Table S{k}'),
        f'Enriched GO terms in the {c} condition (revised: corrected fitness values and current GO annotation).')
add('Table S13', rows_of(R('Table_S13_revised.xlsx'), 'Table S13'),
    'Enriched GO processes in glucose with the largest proteome mass (top 15 per category), each assigned to a single '
    'primary biological process (revised).', header_row=3)
add('Table S13 all terms', rows_of(R('Table_S13_revised.xlsx'), 'All glucose terms'), 'All enriched GO terms in glucose with their mass.')
add('Table S14', rows_of(R('Table_S14_revised.xlsx'), 'Table S14'),
    'GO terms enriched among the 1,946 core mean-effect genes (revised, current GO annotation).', header_row=3)
add('Table S14 core fraction', rows_of(R('Table_S14_revised.xlsx'), 'Core fraction by condition'),
    'Proteome fraction of the core mean-effect genes by condition.')
add('Table S15', rows_of(R('Table_S15_revised.xlsx'), 'Table S15'),
    'GO terms enriched among the 108 recurrent fitness-enhancing candidates (revised: current GO annotation), with core genes and probable Tn-seq '
    'artefacts (rpsG, rpsR) marked (revised).', header_row=2)
add('Table S16', [list(r) for r in orig['Table S16'].iter_rows(values_only=True) if any(v is not None for v in r)],
    'Genes of each sigma-factor regulon (RegulonDB; unchanged).')
add('Table S17', rows_of(R('Table_S17_S18_revised.xlsx'), 'Table S17'),
    'High-confidence function predictions for poorly annotated fitness-enhancing genes (revised: core column).')
add('Table S18', rows_of(R('Table_S17_S18_revised.xlsx'), 'Table S18'),
    'Low-confidence function predictions for poorly annotated fitness-enhancing genes (revised: core column).')
add('Table S19', rows_of(R('Table_S19_revised.xlsx'), 'Table S19'),
    'Number of conditions (of 46) in which each gene falls in each fitness category, with the recurrent fitness-'
    'enhancing candidates, the statistically supported core and probable artefacts (revised).', header_row=2)

# S20: statistical tests (one sheet per test) ---------------------------------------------------------------------
s20 = openpyxl.load_workbook(R('Table_S20_statistical_tests.xlsx'), read_only=True)
for n in s20.sheetnames:
    name = 'S20 index' if n == 'Index' else n
    ws = wb.create_sheet(name[:31])
    for r in s20[n].iter_rows(values_only=True):
        ws.append(list(r))
    ws['A1'].font = BOLD
readme.append(['Table S20', 'Statistical tests added in the revision, one sheet per test (S20a-S20w); see sheet "Data 20 index".'])

# S21-S25: new tables ----------------------------------------------------------------------------------------------
g = pd.read_excel(R('Table_SX_essential_without_fitness.xlsx'), sheet_name='Genes (287)', header=1)
add('Table S21', df_rows(g),
    '287 genes without RB-TnSeq fitness data annotated as essential in PEC, DEG or the Keio collection; added to the '
    'essential category in the GO enrichment only.',
    title='Supplementary Table S21. Genes without fitness data annotated as essential (PEC, Hashimoto et al. 2005; DEG, '
          'Liang et al. 2024; Keio, Yamamoto et al. 2009), used only in the GO enrichment of the essential category.')
t = pd.read_csv(R('R1_analyses', 'output', 'a28_fig_S6_terms_themes.csv'))
add('Table S22', df_rows(t),
    'GO biological processes enriched among mean-effect and fitness-enhancing genes in the nine species, grouped into '
    'functional themes (Supplementary Figure S7B).',
    title='Supplementary Table S22. Enriched GO processes of the multi-species analysis and their functional theme '
          '(topGO elim + Fisher; p < 0.001, gene frequency > 0.3, >= 2 annotated genes).')
pp = pd.read_csv(os.path.join(HERE, 'external', 'Pp_FigureS5_share', 'results', 'FigureS5_tableC.tsv'), sep='\t')
add('Table S23', df_rows(pp),
    'The 168 mean-effect genes among the top 5% most expressed genes of P. putida KT2440 in glucose (Supplementary '
    'Figure S6A).',
    title='Supplementary Table S23. Highly expressed mean-effect genes of Pseudomonas putida KT2440 in glucose '
          '(fitness: Thompson et al. 2019, experiment set6IT057; expression: GSE85475, log10(TPM + 1)).')
cog = pd.read_csv(os.path.join(HERE, 'external', 'Pp_FigureS5_share', 'results', 'TableSX_COG_enrichment.tsv'), sep='\t')
add('Table S24', df_rows(cog),
    'COG enrichment of the highly expressed mean-effect genes of P. putida versus the remaining mean-effect genes '
    '(Supplementary Figure S6B).',
    title="Supplementary Table S24. COG categories of the 168 highly expressed mean-effect genes of P. putida compared "
          "with the remaining mean-effect genes (Fisher's exact test, Benjamini-Hochberg FDR).")
lof = pd.read_csv(R('R1_analyses', 'output', 'a09_all_ltee_lof_genes.csv'))
cand = pd.read_csv(R('R1_analyses', 'output', 'a05_enhancing_108_candidates_32_core.csv'))
lof['enhancing_candidate_108'] = lof.bnumber.isin(cand.bnumber)
lof['enhancing_core_71'] = lof.bnumber.isin(cand.loc[cand.significant_core, 'bnumber'])
lof = lof.drop(columns=['enhancing_core_32'])
add('Table S25', df_rows(lof),
    'Genes with loss-of-function mutations in the 12 LTEE populations (ALEdb), with their glucose fitness category '
    '(Supplementary Figure S5).',
    title='Supplementary Table S25. Genes with loss-of-function mutations in the 12 LTEE populations (Ara+1-6, Ara-1-6). '
          'in_universe: gene with glucose fitness data (433 genes); category: fitness category in glucose; '
          'category_t4: category requiring |t| > 4.')
readme.column_dimensions['A'].width = 22; readme.column_dimensions['B'].width = 120

# npj format (28 Oct 2026): the Excel tables are "Supplementary Data <n>" (same numbers; S20b -> 20b). Sheet names,
# titles and the ReadMe are renamed; figure panels in lower case. Data cells (genes, values) are not touched.
import re
from copy import copy


def rename_text(v):
    v = re.sub(r'Supplementary Tables? S(\d+)', r'Supplementary Data \1', v)
    v = re.sub(r'\bTables? S(\d+)', r'Supplementary Data \1', v)
    v = re.sub(r'^S(20[a-w]\d?)_[A-Za-z0-9_]+\. ', r'Supplementary Data \1. ', v)   # title of each S20 sheet
    v = re.sub(r'\bS(20[a-w]\d?)\b', r'\1', v)                                     # S20b, S20a-S20w in the index
    v = re.sub(r'\bS(20[a-w]\d?_)', r'\1', v)                                       # sheet names listed in the index
    v = re.sub(r'(Figures? S\d+)([A-H])\b', lambda m: m.group(1) + m.group(2).lower(), v)
    return re.sub(r'(Fig\. \d)([A-H])\b', lambda m: m.group(1) + m.group(2).lower(), v)


def rename_sheet(n):
    n = re.sub(r'^Table S(\d+)', r'Data \1', n)
    n = re.sub(r'^S20 index$', 'Data 20 index', n)
    return re.sub(r'^S(20[a-w]\d?)_', r'\1_', n)[:31]


DATA_OF = {}
for ws in wb.worksheets:
    if ws.title in ('ReadMe', 'S20 index'):
        if ws.title == 'S20 index':
            ws.title = 'Data 20 index'; DATA_OF[ws.title] = 20
        for row in ws.iter_rows():
            for c in row:
                if isinstance(c.value, str):
                    c.value = rename_text(c.value)
        continue
    k = int(re.match(r'(?:Table S|S)(\d+)', ws.title).group(1))
    for row in ws.iter_rows(max_row=3):
        for c in row:
            if isinstance(c.value, str):
                c.value = rename_text(c.value)
    ws.title = rename_sheet(ws.title)
    DATA_OF[ws.title] = k
readme['A1'] = 'Supplementary Data'
wb.save(OUT)

# one file per Supplementary Data (npj: "Supplementary Data XX"): values, bold, wrap, widths and frozen panes copied
SD = R('Supplementary_Data'); os.makedirs(SD, exist_ok=True)
for k in sorted(set(DATA_OF.values())):
    w = openpyxl.Workbook(); w.remove(w.active)
    for ws in wb.worksheets:
        if DATA_OF.get(ws.title) != k:
            continue
        ns = w.create_sheet(ws.title)
        for row in ws.iter_rows():
            for c in row:
                if c.value is not None:
                    d = ns.cell(row=c.row, column=c.column, value=c.value)
                    if c.font is not None and c.font.b:
                        d.font = copy(c.font)
                    if c.alignment is not None and c.alignment.wrap_text:
                        d.alignment = copy(c.alignment)
        for col, dim in ws.column_dimensions.items():
            ns.column_dimensions[col].width = dim.width
        ns.freeze_panes = ws.freeze_panes
    w.save(os.path.join(SD, f'Supplementary_Data_{k}.xlsx'))
print('Supplementary Data files:', len(set(DATA_OF.values())))
print('saved', OUT, len(wb.sheetnames), 'sheets;', 'S25 genes in universe:', int(lof.in_universe.sum()), 'of', len(lof))
