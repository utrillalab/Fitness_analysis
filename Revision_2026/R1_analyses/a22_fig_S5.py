"""Revised Supplementary Figure S5 (R1-M7a, M7b): P. putida KT2440, fitness on glucose vs expression.

Same data and selection as O. A. Aguilar-Vera's notebook (../external/Pp_FigureS5_share; reproduced in a21_putida_reproduce.py):
  A  fitness (set6IT057 D-Glucose) vs log10(TPM + 1) (GSE85475 controls), coloured by fitness category, with marginal
     distributions; dashed line = top 5% expression cutoff. The 168 highly expressed mean-effect genes are outlined.
  B  COG composition of the highly expressed mean-effect genes vs the remaining mean-effect genes (single-letter COGs),
     Fisher's exact test per COG letter, BH correction (* FDR < 0.05; fold enrichment shown).
Larger fonts and vector output for M7b. Category colours as Supplementary Figures 1 and 4.
Output: output/Fig_S5_revised.{svg,pdf,png}
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import fisher_exact
from common import *

plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
COL = {'essential': '#b2182b', 'important': '#e08214', 'mean-effect': '#bdbdbd', 'enhancing': '#2166ac'}
LAB = {'essential': 'Essential', 'important': 'Important', 'mean-effect': 'Mean-effect', 'enhancing': 'Enhancing'}
INK, INK2 = '#222222', '#666666'
TOPC, RESTC = '#4d4d4d', '#d9d9d9'

PP = os.path.join(os.path.dirname(HERE), 'external', 'Pp_FigureS5_share')
COG = pd.DataFrame({'COG': list('JAKLBDYVTMNZWUOXCGEFHIPQRS'),
                    'group': ['Information'] * 5 + ['Cellular processes'] * 11 + ['Metabolism'] * 8 + ['Poorly characterized'] * 2,
                    'name': ['Translation', 'RNA processing', 'Transcription', 'Replication and repair', 'Chromatin',
                             'Cell division', 'Nuclear structure', 'Defense mechanisms', 'Signal transduction',
                             'Cell wall/membrane', 'Cell motility', 'Cytoskeleton', 'Extracellular structures',
                             'Intracellular trafficking', 'Posttranslational mod., chaperones', 'Mobilome',
                             'Energy production', 'Carbohydrate metabolism', 'Amino acid metabolism',
                             'Nucleotide metabolism', 'Coenzyme metabolism', 'Lipid metabolism', 'Inorganic ion transport',
                             'Secondary metabolites', 'General function', 'Function unknown']}).set_index('COG')

# ---- data (as a21) ----------------------------------------------------------------------------------
tpm = pd.read_csv(os.path.join(PP, 'data', 'processed', 'GSE85475_TPM.tsv'), sep='\t')
tpm['logTPM'] = np.log10(tpm[['TPM_C_1_0', 'TPM_C_2_0', 'TPM_C_3_0']].mean(axis=1) + 1)
fit = pd.read_csv(os.path.join(PP, 'data', 'raw', 'fit_organism_Putida.tsv'), sep='\t')
d = fit[['locusId', 'set6IT057 D-Glucose (C)']].rename(columns={'set6IT057 D-Glucose (C)': 'fitness'}) \
    .merge(tpm[['locus_tag', 'logTPM']], left_on='locusId', right_on='locus_tag')
d['category'] = classify(d.fitness)
cut = d.logTPM.quantile(0.95)
d['top'] = d.logTPM >= cut
ft = pd.read_table(os.path.join(PP, 'data', 'raw', 'GCA_000007565_feature_table.txt'))
cogs = pd.read_table(os.path.join(PP, 'data', 'raw', 'COGs.csv')).rename(columns={'bnum': 'locus_tag'})
genes = ft.loc[ft['# feature'] == 'gene', ['locus_tag']].drop_duplicates().merge(cogs[['locus_tag', 'func']], how='left')
me = d[d.category == 'mean-effect'].merge(genes, how='left', on='locus_tag')
me1 = me[me.func.notna() & ~me.func.str.contains(',', na=False)]
nt, nr = int(me1.top.sum()), int((~me1.top).sum())
rows = []
for k in me1.func.unique():
    a = int(((me1.func == k) & me1.top).sum()); b = int(((me1.func == k) & ~me1.top).sum())
    o, p = fisher_exact([[a, nt - a], [b, nr - b]])
    rows.append({'COG': k, 'top_n': a, 'rest_n': b, 'top_%': 100 * a / nt, 'rest_%': 100 * b / nr, 'odds': o, 'p': p})
e = pd.DataFrame(rows).set_index('COG').join(COG)
e['FDR'] = bh(e.p.values)
e['fold'] = e['top_%'] / e['rest_%']
GORDER = ['Information', 'Cellular processes', 'Metabolism', 'Poorly characterized']
e['g'] = e.group.map({g: i for i, g in enumerate(GORDER)})
e = e.sort_values(['g', 'rest_%'], ascending=[True, False])
gsum = e.groupby('group')[['top_n', 'rest_n']].sum()

# ---- figure ------------------------------------------------------------------------------------------
fig = plt.figure(figsize=(7.2, 8.4))
gtop = GridSpec(2, 2, figure=fig, height_ratios=[0.22, 1], width_ratios=[1, 0.16], hspace=0.04, wspace=0.03,
                left=0.10, right=0.97, top=0.95, bottom=0.60)
gbot = GridSpec(1, 1, figure=fig, left=0.36, right=0.95, top=0.49, bottom=0.06)
ax = fig.add_subplot(gtop[1, 0]); axx = fig.add_subplot(gtop[0, 0], sharex=ax); axy = fig.add_subplot(gtop[1, 1], sharey=ax)

XMIN, XMAX = -7, 2
for k in ['mean-effect', 'important', 'essential', 'enhancing']:
    s = d[d.category == k]
    ax.scatter(s.fitness.clip(XMIN, XMAX), s.logTPM, s=6, c=COL[k], lw=0, alpha=0.8 if k == 'mean-effect' else 0.95,
               label=f'{LAB[k]} (n = {len(s):,})', zorder=2 if k == 'mean-effect' else 3, rasterized=True)
hi = d[(d.category == 'mean-effect') & d.top]
ax.scatter(hi.fitness, hi.logTPM, s=7, facecolor='#636363', edgecolor='black', lw=0.3, zorder=4,
           label=f'Mean-effect, top 5% expressed (n = {len(hi)})', rasterized=True)
for x in (-2, -0.5, 0.5):
    ax.axvline(x, ls=':', c='grey', lw=0.7, zorder=1)
ax.axhline(cut, ls='--', c=INK, lw=0.7, zorder=1)
axy.text(0.95, cut + 0.08, f'Top 5%\n(TPM ≥ {10 ** cut - 1:.0f};\n{int(d.top.sum())} genes)', transform=axy.get_yaxis_transform(),
         fontsize=6.5, color=INK, va='bottom', ha='right')
ax.set_xlim(XMIN, XMAX); ax.set_ylim(-0.1, d.logTPM.max() + 0.2)
ax.set_xlabel('Fitness in glucose (set6IT057; values below −7 shown at −7)'); ax.set_ylabel('Expression, log$_{10}$(TPM + 1)')
ax.legend(frameon=False, fontsize=7, loc='upper left', markerscale=2, handletextpad=0.2, borderaxespad=0.1)
for t, x in [('Essential', -4.5), ('Important', -1.25), ('Mean-\neffect', 0), ('Enh.', 1.25)]:
    axx.text(x, 1.0, t, transform=axx.get_xaxis_transform(), ha='center', va='bottom', fontsize=7, color=INK2)
axx.hist(d.fitness.clip(XMIN, XMAX), bins=np.arange(XMIN, XMAX + 0.05, 0.05), color='#969696')
axx.set_yscale('log'); axx.axis('off')
axy.hist(d.logTPM, bins=60, orientation='horizontal', color='#969696')
axy.axhline(cut, ls='--', c=INK, lw=0.7); axy.axis('off')
ax.text(-0.11, 1.25, 'a', transform=ax.transAxes, fontsize=11, fontweight='bold')

axb = fig.add_subplot(gbot[0, 0])
y = np.arange(len(e))[::-1].astype(float)
# gap between groups
off, ys, prev = 0, [], None
for g in e.group:
    if prev is not None and g != prev:
        off += 0.8
    ys.append(off); off += 1; prev = g
y = -np.array(ys)
h = 0.38
axb.barh(y + h / 2, e['top_%'], h, color=TOPC, label=f'Mean-effect, top 5% expressed (n = {nt})')
axb.barh(y - h / 2, e['rest_%'], h, color=RESTC, label=f'Other mean-effect genes (n = {nr:,})')
for yi, (_, r) in zip(y, e.iterrows()):
    if r.FDR < 0.05:
        axb.text(max(r['top_%'], r['rest_%']) + 0.3, yi, f'{r.fold:.1f}× *', va='center', fontsize=7.5, fontweight='bold', color=INK)
axb.set_yticks(y, [f'{n} ({c})' for c, n in zip(e.index, e.name)], fontsize=7.5)
axb.tick_params(axis='y', length=0)
axb.set_xlabel('% of genes with a single COG category')
axb.legend(frameon=False, fontsize=7.5, loc='lower right')
xmax = e[['top_%', 'rest_%']].max().max() * 1.25
axb.set_xlim(0, xmax)
prev = None
for yi, g in zip(y, e.group):
    if g != prev:
        gy = y[e.group.values == g]
        tn, rn = gsum.loc[g, 'top_n'], gsum.loc[g, 'rest_n']
        gl = 'Poorly char.' if g == 'Poorly characterized' else g
        axb.text(-0.53, gy.mean(), f'{gl}\n{100 * tn / nt:.0f}% vs {100 * rn / nr:.0f}%', transform=axb.get_yaxis_transform(),
                 ha='center', va='center', fontsize=7.5, fontweight='bold', color=INK, rotation=90)
        axb.plot([-0.49, -0.49], [gy.min() - 0.3, gy.max() + 0.3], transform=axb.get_yaxis_transform(), c=INK2, lw=0.6,
                 clip_on=False)
        prev = g
axb.text(-0.6, 1.04, 'b', transform=axb.transAxes, fontsize=11, fontweight='bold')
axb.set_title('COG categories of highly expressed mean-effect genes (* FDR < 0.05, Fisher\'s exact test)',
              fontsize=8, loc='left', x=-0.53, pad=8)

for ext in ('svg', 'pdf'):
    fig.savefig(f'{OUT}/Fig_S5_revised.{ext}', dpi=600)
fig.savefig(f'{OUT}/Fig_S5_revised.png', dpi=300)
print(e[['name', 'top_n', 'rest_n', 'fold', 'FDR']].round(3).to_string())
print('group totals (top % vs rest %):', {g: (round(100 * gsum.loc[g, 'top_n'] / nt), round(100 * gsum.loc[g, 'rest_n'] / nr)) for g in GORDER})
