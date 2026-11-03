"""Regenerate the fitness input files of utrillalab/Fitness_analysis with correctly paired replicates.

Problem: allgen4691_fitfg_newname.csv (= Supplementary Table S1), genes_categories/genesfitness_*.csv
and all_genes_categorized.csv average replicate pairs from different conditions (off-by-one pairing).
Fix: take fitness as the mean of the two replicates of the SAME condition, exactly as
`fitness_averages` in Code_Notebooks/Final_Table_and_Proteome_Fractions.ipynb (cell 13).

Usage:  python regenerate_inputs.py <path to Fitness_analysis clone> <output dir>
"""
import sys, os
import numpy as np
import pandas as pd

repo = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'regenerated_inputs')
os.makedirs(os.path.join(out, 'genes_categories'), exist_ok=True)
inf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'original_submission', 'Infiles')  # files as first submitted

C = ['Glucose', 'LB', 'Glucosamine', 'Glycerol', 'Xylose', 'Succinate', 'Mannose', 'Galactose',
     'Fructose', 'Acetate', 'Pyruvate']

# --- correct fitness averages (copy of notebook cell 13) ---------------------------------------
ft = pd.read_csv(os.path.join(repo, 'Input_Data', 'fit_organism_Keio.tsv'), sep='\t')
fa = pd.DataFrame(index=ft.sysName)
for c in range(5, 95, 2):
    assert ft.columns[c].split(' ', 1)[1] == ft.columns[c + 1].split(' ', 1)[1], 'replicates differ'
    fa[ft.columns[c]] = ft.iloc[:, [c, c + 1]].mean(axis=1).values
fa['LB'] = ft.iloc[:, [112, 131]].mean(axis=1).values
col = {'Glucose': 'set1IT003 D-Glucose (C)', 'Acetate': 'set1IT027 acetate (C)',
       'Glucosamine': 'set1IT063 D-Glucosamine Hydrochloride (C)', 'Glycerol': 'set1IT047 Glycerol (C)',
       'Pyruvate': 'set1IT037 pyruvate (C)', 'Xylose': 'set1IT011 D-Xylose (C)',
       'Mannose': 'set1IT065 D-Mannose (C)', 'Galactose': 'set1IT013 D-Galactose (C)',
       'Succinate': 'set1IT039 succinate (C)', 'Fructose': 'set1IT005 D-Fructose (C)', 'LB': 'LB'}
fit = pd.DataFrame({c: fa[col[c]] for c in C})
# Some repo files use older b-numbers (e.g. ybbD b0501 vs b0500 in the raw table): fall back to gene name.
fit_by_name = fit.copy()
fit_by_name.index = ft.set_index('sysName').loc[fit.index, 'geneName'].values
fit_by_name = fit_by_name[~fit_by_name.index.duplicated()]


def lookup(df, c):
    by_b = df['bnum'].map(fit[c])
    return by_b.fillna(df['gene'].map(fit_by_name[c]))


def categ(f):
    """Manuscript / Final_Table notebook rule: mean-effect is -0.5 <= f <= 0.5 (inclusive)."""
    if pd.isna(f):
        return np.nan
    if f <= -2:
        return 'Essen'
    if -0.5 <= f <= 0.5:
        return 'No_effect'
    if f > 0.5:
        return 'Enhancing'
    return 'Important'


long_label = {'Essen': 'Essential', 'Important': 'Important', 'No_effect': 'No effect', 'Enhancing': 'Enhancing'}

# --- 1. allgen4691_fitfg_newname.csv (also the basis of Table S1) ------------------------------
ag = pd.read_csv(os.path.join(inf, 'allgen4691_fitfg_newname.csv'))
old_ag = ag.copy()
for c in C:
    ag[c + '_fit'] = lookup(ag, c)          # genes absent from RB-TnSeq -> NaN (was 0)
ag.to_csv(os.path.join(out, 'allgen4691_fitfg_newname.csv'), index=False, na_rep='NA')

# --- 2. genes_categories/genesfitness_<cond>.csv  and  3. all_genes_categorized.csv ------------
report, changed_rows, long_rows = [], [], []
for c in C:
    old = pd.read_csv(os.path.join(inf, 'genes_categories', f'genesfitness_{c}.csv'))
    new = old[['bnum', 'gene']].copy()
    new['fitval'] = lookup(new, c)
    new['Categ'] = new['fitval'].map(categ)
    missing = new['fitval'].isna().sum()
    new = new.dropna(subset=['fitval'])
    new.to_csv(os.path.join(out, 'genes_categories', f'genesfitness_{c}.csv'), index=False)

    m = old.merge(new, on=['bnum', 'gene'], suffixes=('_old', '_new'))
    diff = m[m.Categ_old != m.Categ_new]
    # how many of the changes are only the +0.5 boundary convention of the old files
    boundary = ((m.fitval_old == 0.5) & (m.Categ_old == 'Enhancing')).sum()
    counts_old = old.Categ.value_counts()
    counts_new = new.Categ.value_counts()
    report.append({'condition': c, 'genes': len(new), 'changed_category': len(diff),
                   'changed_%': round(100 * len(diff) / len(new), 1), 'old_boundary_0.5': int(boundary),
                   **{f'{k}_old': int(counts_old.get(k, 0)) for k in long_label},
                   **{f'{k}_new': int(counts_new.get(k, 0)) for k in long_label},
                   'genes_missing_in_raw': int(missing)})
    d = diff[['bnum', 'gene', 'fitval_old', 'Categ_old', 'fitval_new', 'Categ_new']].copy()
    d.insert(0, 'condition', c)
    changed_rows.append(d)
    lr = new.rename(columns={'fitval': 'value'})
    lr['condition'] = c
    lr['category'] = lr['Categ'].map(long_label)
    long_rows.append(lr[['bnum', 'gene', 'condition', 'category', 'value']])

pd.concat(long_rows).to_csv(os.path.join(out, 'all_genes_categorized.csv'), index=False)
rep = pd.DataFrame(report).set_index('condition')
rep.to_csv(os.path.join(out, 'change_report_by_condition.csv'))
pd.concat(changed_rows).to_csv(os.path.join(out, 'genes_changing_category.csv'), index=False)

# --- 4. Corrected Supplementary Table S1 (same layout: fitness block then proteomics block) ----
s1 = ag[['bnum', 'gene'] + [c + '_fit' for c in C] + C].copy()
s1.columns = ['bnum', 'gene'] + [f'Fitness_{c}' for c in C] + [f'Proteome_fg_{c}' for c in C]
s1.to_csv(os.path.join(out, 'Table_S1_corrected.csv'), index=False, na_rep='NA')

print(rep[['genes', 'changed_category', 'changed_%', 'old_boundary_0.5', 'genes_missing_in_raw']].to_string())
print('\nCategory counts old -> new:')
for k in long_label:
    print(f'  {k:10s}', ' '.join(f'{c[:4]} {r[k + "_old"]}->{r[k + "_new"]}' for c, r in rep.iterrows()))
print('\nTable S1 genes with NA fitness:', int(ag[C[0] + '_fit'].isna().sum()),
      '| previously coded 0 in all conditions:', int((old_ag[[c + '_fit' for c in C]] == 0).all(axis=1).sum()))
print('Output written to', out)
