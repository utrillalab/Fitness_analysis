"""R1-m16: summary panel of the GO processes enriched in mean-effect and fitness-enhancing genes across nine species.

Input:
  - ../external/Liz_GO_nine_species/all_gos.csv: topGO results of Liz's notebook EnrichmentGO_several_orgs_Suppl_S6.ipynb (glucose; elim + Fisher;
    UniProt GO; top 20 terms per species and category), exported on 14 Oct 2026. B. thetaiotaomicron is excluded (not in
    the text).
  - E. coli: corrected ViSEAGO glucose results (output/viseago/corrected/BP_Results_Glucose.xls).
One filter for every species (as in Fig. 3): p < 0.001, gene frequency > 0.3, >= 2 annotated genes. Each term is assigned
to a functional theme by keywords (first match in THEMES); the cells give the number of enriched terms per theme.
Output: output/Fig_S6_summary.png/svg/pdf and output/a28_fig_S6_terms_themes.csv (term -> theme, for a supplementary table).
"""
import os, re
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import OUT, HERE

P, FREQ, NANN = 0.001, 0.3, 2
SPECIES = {'escherichia_coli': 'E. coli', 'pseudomonas_putida_kt2440': 'P. putida',
           'pseudomonas_fluorescens_fw300-n2c3': 'P. fluorescens', 'pseudomonas_simiae_wcs417': 'P. simiae',
           'pseudomonas_stutzeri_rch2': 'P. stutzeri', 'caulobacter_crescentus_na1000': 'C. crescentus',
           'dinoroseobacter_shibae_dfl_12': 'D. shibae', 'sinorhizobium_meliloti_1021': 'S. meliloti',
           'echinicola_vietnamensis': 'E. vietnamensis'}
THEMES = [  # (theme, keywords); first match wins
    ('Motility, chemotaxis and adhesion', r'flagell|chemotaxis|motility|taxis|pilus|adhesion|biofilm'),
    ('Proteolysis and protein secretion', r'proteolysis|peptid(ase|e catabolic)|protein secretion|protein localization|protein transport'),
    ('Transport', r'transport|export|import|secretion'),
    ('Catabolism', r'catabolic|degradation'),
    ('Biosynthesis', r'biosynthe'),
    ('Signalling and regulation', r'signal|phosphorelay|regulation|transcription|communication'),
    ('Stress response and DNA repair', r'stress|DNA repair|DNA damage|response to|SOS'),
    ('Other metabolic processes', r'metabolic|metabolism'),
]
CATS = [('noeffect', 'Mean-effect'), ('enhance', 'Fitness-enhancing')]


def theme(term):
    for name, pat in THEMES:
        if re.search(pat, term, re.I):
            return name
    return 'Other'


# other species (Liz's topGO results)
g = pd.read_csv(os.path.join(os.path.dirname(HERE), 'external', 'Liz_GO_nine_species', 'all_gos.csv'))
g = g[g.Organism.isin(SPECIES) & g.Category.isin([c for c, _ in CATS])]
g = g.assign(freq=g.Significant / g.Annotated)[['Organism', 'Category', 'GO', 'Term', 'Annotated', 'Significant', 'freq', 'PValue']]

# E. coli (corrected ViSEAGO, glucose)
r = pd.read_csv(os.path.join(OUT, 'viseago', 'corrected', 'BP_Results_Glucose.xls'), sep='\t')
rows = []
for col, cat in [('No_effect', 'noeffect'), ('Fitness_enhancing', 'enhance')]:
    for _, x in r[r[f'{col}.pvalue'].notna()].iterrows():
        sig, ann = map(int, re.search(r'\((\d+)/(\d+)\)', x[f'{col}.genes_frequency']).groups())
        rows.append({'Organism': 'escherichia_coli', 'Category': cat, 'GO': x['GO.ID'], 'Term': x['term'],
                     'Annotated': ann, 'Significant': sig, 'freq': sig / ann if ann else 0, 'PValue': x[f'{col}.pvalue']})
d = pd.concat([g, pd.DataFrame(rows)], ignore_index=True)
d = d[(d.PValue < P) & (d.freq > FREQ) & (d.Annotated >= NANN)].copy()
d['Theme'] = d.Term.map(theme)
d['Species'] = d.Organism.map(SPECIES)
d['Category'] = d.Category.map(dict(CATS))
d.sort_values(['Category', 'Species', 'Theme', 'PValue'])[
    ['Category', 'Species', 'GO', 'Term', 'Theme', 'Annotated', 'Significant', 'PValue']].to_csv(
    f'{OUT}/a28_fig_S6_terms_themes.csv', index=False)

themes = [t for t, _ in THEMES] + ['Other']
sp = list(SPECIES.values())
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True, gridspec_kw={'wspace': 0.08})
vmax = max(d.groupby(['Category', 'Species', 'Theme']).size().max(), 1)
for ax, (_, cname) in zip(axes, CATS):
    m = (d[d.Category == cname].groupby(['Theme', 'Species']).size().unstack(fill_value=0)
         .reindex(index=themes, columns=sp, fill_value=0))
    ax.imshow(np.ma.masked_equal(m.values, 0), cmap='Blues', vmin=0, vmax=vmax * 1.15, aspect='auto')
    ax.set_facecolor('#f4f4f4')
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.values[i, j]
            if v:
                ax.text(j, i, v, ha='center', va='center', fontsize=9,
                        color='white' if v > vmax * 0.55 else '#1a1a1a')
    ax.set_xticks(range(len(sp))); ax.set_xticklabels(sp, rotation=45, ha='right', fontsize=9, style='italic')
    ax.set_yticks(range(len(themes))); ax.set_yticklabels(themes, fontsize=9)
    ax.set_xticks(np.arange(-.5, len(sp)), minor=True); ax.set_yticks(np.arange(-.5, len(themes)), minor=True)
    ax.grid(which='minor', color='white', lw=2); ax.tick_params(which='both', length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    n_sp = (m.sum(axis=0) > 0).sum()
    ax.set_title(f'{cname} genes ({n_sp} of {len(sp)} species with enriched terms)', fontsize=10)
fig.savefig(f'{OUT}/Fig_S6_summary.png', dpi=300, bbox_inches='tight')
for ext in ('svg', 'pdf'):
    fig.savefig(f'{OUT}/Fig_S6_summary.{ext}', bbox_inches='tight')
print(d.groupby(['Category', 'Theme']).agg(terms=('Term', 'size'), species=('Species', 'nunique')).to_string())
print('\nOther:', d[d.Theme == 'Other'][['Species', 'Category', 'Term']].to_string(index=False))
