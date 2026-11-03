"""Redrawn Fig. 4 (R1-M5a, M5b, m7): enriched GO processes with the largest proteome mass, glucose and galactose.

Changes from the original figure:
  - corrected GO enrichment (output/viseago/corrected/, fitness with correctly paired replicates);
  - mass of a process = summed mass of its genes in that category / total measured proteome
    (no division by the number of GO annotations, R1-m7);
  - redundant terms removed: terms are taken in decreasing order of mass and a term is skipped if it shares
    more than half of its genes (overlap coefficient > 0.5) with a term already shown;
  - top 5 per category (top 2 for fitness-enhancing), as in the original."""
import os, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *

VIS = os.path.join(OUT, 'viseago', 'corrected')
CAT_COL = {'Essential': 'essential', 'Important': 'important', 'No_effect': 'mean-effect',
           'Fitness_enhancing': 'enhancing'}
LABEL = {'essential': 'Essential', 'important': 'Important', 'mean-effect': 'Mean-effect', 'enhancing': 'Enhancing'}
PAL = {'essential': '#b2182b', 'important': '#e08214', 'mean-effect': '#969696', 'enhancing': '#2166ac'}
TOP = {'essential': 5, 'important': 5, 'mean-effect': 5, 'enhancing': 2}

fit, names = load_fitness()
mass = load_mass()
eco = pd.read_csv(os.path.join(REPO, 'Input_Data', 'All-genes-of-E.-coli-K-12-substr-MG1655.txt'), sep='\t')
sym2b = {}
for n, b in zip(eco['Gene Name'], eco['Accession-1']):
    if isinstance(n, str) and isinstance(b, str):
        sym2b.setdefault(n.lower(), b)
for b, n in names.items():
    if isinstance(n, str):
        sym2b.setdefault(n.lower(), b)


def terms_for(cond):
    res = pd.read_csv(os.path.join(VIS, f'BP_Results_{cond}.xls'), sep='\t')
    m = mass[cond].dropna(); tot = m.sum()
    rows = []
    for col, cat in CAT_COL.items():
        sig = res[res[f'{col}.Significant_genes_symbol'].notna()]
        for _, r in sig.iterrows():
            bs = {sym2b.get(s.lower()) for s in str(r[f'{col}.Significant_genes_symbol']).split(';') if s}
            bs.discard(None)
            rows.append({'condition': cond, 'category': cat, 'GO.ID': r['GO.ID'], 'term': r['term'],
                         'n_genes': len(bs), 'n_with_proteomics': int(m.index.isin(list(bs)).sum()), 'genes': bs,
                         'mass_%': 100 * m.reindex(list(bs)).dropna().sum() / tot,
                         '-log10p': r[f'{col}.-log10_pvalue']})
    return pd.DataFrame(rows)


def select(t, n):
    chosen = []
    for _, r in t.sort_values('mass_%', ascending=False).iterrows():
        if any(len(r.genes & c.genes) / max(1, min(len(r.genes), len(c.genes))) > 0.5 for c in chosen):
            continue
        chosen.append(r)
        if len(chosen) == n:
            break
    return pd.DataFrame(chosen)


panels = {}
for c in ['Glucose', 'Galactose']:
    t = terms_for(c)
    panels[c] = pd.concat([select(t[t.category == k], TOP[k]) for k in CATS if (t.category == k).any()])
out = pd.concat(panels.values()).drop(columns='genes')
out.to_csv(f'{OUT}/a20_fig4_terms.csv', index=False)
pd.set_option('display.width', 200)
print(out[['condition', 'category', 'term', 'n_genes', 'n_with_proteomics', 'mass_%']].round(2).to_string(index=False))

# npj style (28 Oct 2026): Arial 7 pt at print size (140 mm wide), colours of npj_style.py (the same as before)
sys.path.insert(0, os.path.dirname(HERE)); from npj_style import *
fig, axes = plt.subplots(2, 1, figsize=(140 * MM, 170 * MM))
xmax = max(p['mass_%'].max() for p in panels.values()) * 1.12
for ax, (c, p), lt in zip(axes, panels.items(), 'ab'):  # npj: lower-case bold panel letters
    p = p.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(p))
    ax.barh(y, p['mass_%'], color=[PAL[k] for k in p.category], edgecolor='black', lw=0.4)
    ax.set_yticks(y); ax.set_yticklabels([f"{t} [{n}]" for t, n in zip(p.term, p.n_with_proteomics)])
    for i in range(1, len(p)):
        if p.category[i] != p.category[i - 1]:
            ax.axhline(i - 0.5, color='grey', ls='--', lw=0.5)
    ax.set_xlim(0, xmax); ax.set_xticks(range(0, int(xmax) + 1, 5)); ax.set_xlabel('Proteome mass fraction (%)')
    ax.set_ylim(-0.6, len(p) - 0.4)
    ax.set_title(c)
handles = [plt.Rectangle((0, 0), 1, 1, color=PAL[k], ec='black', lw=0.4) for k in CATS]
fig.legend(handles, [LABEL[k] for k in CATS], loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0),
           handlelength=1, handleheight=1)
fig.tight_layout(rect=(0, 0, 1, 0.965), h_pad=2)
for ax, lt in zip(axes, 'ab'):
    letter(fig, ax, lt, dx=0)
save(fig, f'{OUT}/Fig04_revised')
