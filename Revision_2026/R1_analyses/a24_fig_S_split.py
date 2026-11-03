"""Revised supplementary figures, split in two (decision of 11 Oct 2026; final numbers to be set when the supplement is
renumbered by order of citation):

Fig_S_classification (fitness classification; provisional S1):
  A  fitness distribution in each condition and pooled, with mode, mean and cut-offs (R1-M2)
  B  fitness vs combined t statistic in glucose (Price et al. 2018 data) (R2-M2a)
  C  mean-effect proteome fraction with cut-offs only and with |t| > 4 (R2-M2a)
  D  proteome mass by category including measured proteins without fitness data (R1-M1)
Fig_S_mass (protein mass and fitness categories; provisional S2, replaces the original S1):
  A  log10 protein mass by category in each condition, Cliff's delta vs mean-effect (original S1, corrected; R1-M4b)
  B  mass distribution by category with quintile bands; C  representation per mass quintile (R1-M4c)
  D  gene-set overlap between conditions; E  conditions per gene (R1-M8)
Panels of the mass figure reuse the drawing functions of a18 and a23 (importing them also refreshes their outputs).
Output: output/Fig_S_classification.{svg,pdf,png}, output/Fig_S_mass.{svg,pdf,png}
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from common import *
import a18_fig_S1_mass_quintiles as q
import a23_fig_S1_revised as s1

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
COL = s1.COL; LAB = s1.LAB; INK, INK2 = s1.INK, s1.INK2
MODEC, MEANC = '#2166ac', '#b2182b'


def letter(ax, s, x=-0.12, y=1.04):
    ax.text(x, y, s.lower(), transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')  # npj: lower-case bold


def save(fig, name):
    for ext in ('svg', 'pdf'):
        fig.savefig(f'{OUT}/{name}.{ext}', dpi=600)
    fig.savefig(f'{OUT}/{name}.png', dpi=300)


# =============================== Fig_S_classification ===============================================
fit, _ = load_fitness()
summ = pd.read_csv(f'{OUT}/a01_distribution_summary.csv').set_index('condition')
P18 = os.path.join(os.path.dirname(HERE), 'Input_Data_Price2018')
tt = pd.read_csv(os.path.join(P18, 'fit_t.tab'), sep='\t').set_index('sysName')
tg = (tt[['set1IT003 D-Glucose (C)', 'set1IT004 D-Glucose (C)']].sum(axis=1) / np.sqrt(2) / np.sqrt(1.5)).reindex(fit.index)
frt = pd.read_csv(f'{OUT}/a05_proteome_fractions_t.csv', header=[0, 1], index_col=0)
miss = pd.read_csv(f'{OUT}/a13_missing_fitness_fractions.csv').set_index('condition')
CORD = ['LB', 'Glycerol', 'Galactose', 'Mannose', 'Glucosamine', 'Pyruvate', 'Fructose', 'Glucose', 'Xylose', 'Succinate',
        'Acetate']
BINS = np.arange(-4, 2.05, 0.05)

fig = plt.figure(figsize=(7.2, 9.8))
gA = GridSpec(3, 4, figure=fig, left=0.09, right=0.98, top=0.96, bottom=0.555, hspace=0.45, wspace=0.18)
for i, c in enumerate(CORD + ['All 11 conditions']):
    ax = fig.add_subplot(gA[i // 4, i % 4])
    v = fit[CORD].values.ravel() if i == 11 else fit[c].values
    v = v[~np.isnan(v)]
    ax.axvspan(-0.5, 0.5, color='#ececec', lw=0, zorder=0)
    ax.hist(np.clip(v, -4, 2), BINS, color='#4d4d4d' if i == 11 else '#969696', lw=0)
    ax.set_yscale('log'); ax.set_xlim(-4.05, 2.05)
    key = c if i < 11 else None
    mode = summ.loc[key, 'mode'] if key else 0.0
    mean = summ.loc[key, 'mean'] if key else np.nanmean(v)
    ax.axvline(mode, color=MODEC, lw=0.9); ax.axvline(mean, color=MEANC, lw=0.9)
    ax.axvline(-2, color=INK, lw=0.6, ls=':')
    me = summ.loc[key, '%mean-effect'] if key else 100 * ((v > -0.5) & (v <= 0.5)).mean()
    ax.set_title(c, fontsize=7.5, pad=2)
    ax.text(0.03, 0.97, f'mean-effect {me:.0f}%\nmean {mean:.2f}', transform=ax.transAxes, fontsize=6, va='top', color=INK)
    ax.set_xticks([-4, -2, 0, 2], ['−4', '−2', '0', '2'] if i >= 8 else [])
    ax.tick_params(axis='y', labelsize=6)
    if i % 4 == 0:
        ax.set_ylabel('Genes')
    if i == 0:
        letter(ax, 'A', -0.32, 1.1)
fig.text(0.53, 0.52, 'Fitness (log$_2$ ratio; values beyond −4 or 2 shown at the edges)', ha='center', fontsize=7.5)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
fig.legend(handles=[Line2D([], [], color=MODEC, lw=1, label='Mode'), Line2D([], [], color=MEANC, lw=1, label='Mean'),
                    Patch(color='#ececec', label='Mean-effect range (−0.5 to 0.5)'),
                    Line2D([], [], color=INK, lw=0.8, ls=':', label='Essential cut-off (−2)')],
           loc='upper center', bbox_to_anchor=(0.53, 0.515), ncol=4, frameon=False, fontsize=6.8, handlelength=1.5)

gB = GridSpec(1, 2, figure=fig, left=0.12, right=0.98, top=0.44, bottom=0.30, wspace=0.32, width_ratios=[1, 1.3])
ax = fig.add_subplot(gB[0, 0])
f, t = fit['Glucose'], tg
sig = t.abs() > 4
ax.scatter(f[~sig].clip(-4, 2), t[~sig].clip(-25, 15), s=3, c='#bdbdbd', lw=0, label='|t| ≤ 4', rasterized=True)
ax.scatter(f[sig].clip(-4, 2), t[sig].clip(-25, 15), s=3, c=MODEC, lw=0, label='|t| > 4', rasterized=True)
for x in (-2, -0.5, 0.5):
    ax.axvline(x, color=INK, lw=0.6, ls=':')
for y in (-4, 4):
    ax.axhline(y, color=MODEC, lw=0.6)
ax.set_xlim(-4.1, 2.1); ax.set_xlabel('Fitness in glucose'); ax.set_ylabel('Combined t statistic')
ax.legend(frameon=False, fontsize=6.5, loc='upper left', markerscale=2.5, handletextpad=0.1)
letter(ax, 'B', -0.3, 1.05)
ax = fig.add_subplot(gB[0, 1])
x = np.arange(len(CORD))
a_ = frt[('cut-offs only', 'mean-effect')].reindex(CORD); b_ = frt[('cut-offs + |t| > 4', 'mean-effect')].reindex(CORD)
ax.bar(x - 0.2, a_, 0.4, color='#636363', label='Cut-offs only')
ax.bar(x + 0.2, b_, 0.4, color=MODEC, label='Cut-offs and |t| > 4')
ax.set_xticks(x, CORD, rotation=45, ha='right'); ax.set_ylim(0, 100)
ax.set_ylabel('Mean-effect proteome\nfraction (%)'); ax.legend(frameon=False, fontsize=6.5, loc='upper right', ncol=2)
letter(ax, 'C', -0.17, 1.05)

gD = GridSpec(1, 1, figure=fig, left=0.09, right=0.8, top=0.21, bottom=0.075)
ax = fig.add_subplot(gD[0, 0])
order = ['mean-effect', 'enhancing', 'important', 'essential']
bottom = np.zeros(len(CORD))
for k in order:
    v = miss[k].reindex(CORD).values
    ax.bar(x, v, 0.7, bottom=bottom, color=COL[k], label=LAB[k], edgecolor='white', lw=0.3)
    if k == 'mean-effect':
        for xi, vi in zip(x, v):
            ax.text(xi, vi / 2, f'{vi:.0f}', ha='center', va='center', fontsize=6, color=INK)
    bottom += v
v = miss['no fitness data'].reindex(CORD).values
ax.bar(x, v, 0.7, bottom=bottom, facecolor='white', edgecolor=COL['essential'], hatch='/////', lw=0.6,
       label='No fitness data')
ax.set_xticks(x, CORD, rotation=45, ha='right'); ax.set_ylim(0, 100); ax.set_ylabel('Measured proteome (%)')
ax.legend(frameon=False, fontsize=6.5, loc='center left', bbox_to_anchor=(1.01, 0.5), reverse=True)
letter(ax, 'D', -0.09, 1.04)
save(fig, 'Fig_S_classification')

# =============================== Fig_S_mass ===========================================================
fig = plt.figure(figsize=(7.2, 9.2))
gA = GridSpec(2, 6, figure=fig, left=0.08, right=0.99, top=0.96, bottom=0.57, hspace=0.55, wspace=0.12)
axs = [fig.add_subplot(gA[i // 6, i % 6]) for i in range(12)]
for a in axs[1:]:
    a.sharey(axs[0])
s1.boxes(axs[:11])
axs[11].axis('off')
axs[11].text(-0.1, 0.5, 'ME, mean-effect\nImp, important\nEss, essential\nEnh, enhancing\n\nNumbers: Cliff\'s delta\n'
             'vs mean-effect\n* Holm-adjusted p < 0.05\n** < 0.01, *** < 0.001\n(Kruskal-Wallis)', fontsize=6.3,
             va='center', color=INK)
for a in axs:
    if a not in (axs[0], axs[6]):
        a.tick_params(labelleft=False)
axs[0].set_ylabel('log$_{10}$ protein mass (fg)'); axs[6].set_ylabel('log$_{10}$ protein mass (fg)')
letter(axs[0], 'A', -0.45, 1.12)

gB = GridSpec(1, 2, figure=fig, left=0.08, right=0.92, top=0.47, bottom=0.27, wspace=0.45, width_ratios=[1, 1.25])
a = fig.add_subplot(gB[0, 0]); b = fig.add_subplot(gB[0, 1])
q.panel_A(a); q.panel_B(b, fig)
letter(a, 'B', -0.2, 1.06); letter(b, 'C', -0.3, 1.06)

gD = GridSpec(1, 3, figure=fig, left=0.115, right=0.98, top=0.19, bottom=0.07, wspace=0.55, width_ratios=[0.75, 1.1, 0.95])
a = fig.add_subplot(gD[0, 0]); b = fig.add_subplot(gD[0, 1]); c = fig.add_subplot(gD[0, 2])
s1.panel_jaccard(a); s1.panel_nconds(b)
b.legend(frameon=False, fontsize=5.8, loc='upper center', handlelength=1.0)
letter(a, 'D', -0.5, 1.04); letter(b, 'E', -0.3, 1.04)

# F: mass of each enriched GO term vs random gene sets of the same size (R1-M5b; a17_go_mass_null.py)
gm = pd.read_csv(f'{OUT}/a17_go_mass_null.csv')
gsum = pd.read_csv(f'{OUT}/a17_go_mass_null_summary.csv').set_index('category')
for k in ['mean-effect', 'important', 'enhancing', 'essential']:
    g = gm[gm.category == k]
    c.scatter(g['null_all_mean_%'], g['observed_mass_%'], s=4, lw=0, alpha=0.7, color=COL[k],
              label=f"{LAB[k]} ({gsum.loc[k, 'median obs/exp vs all']:.1f}×)")
lim = [1e-3, 60]
c.plot(lim, lim, ls='--', color=INK, lw=0.6)
c.set_xscale('log'); c.set_yscale('log'); c.set_xlim(lim); c.set_ylim(lim)
c.set_xlabel('Expected mass, random genes\nof the same number (%)'); c.set_ylabel('Observed mass of the GO term (%)')
c.text(0.03, 0.97, 'Median obs/exp', transform=c.transAxes, fontsize=5.8, va='top', color=INK)
for n, k in enumerate(['essential', 'enhancing', 'important', 'mean-effect']):
    c.text(0.03, 0.90 - 0.07 * n, f"{LAB[k]} {gsum.loc[k, 'median obs/exp vs all']:.1f}×", transform=c.transAxes,
           fontsize=5.8, va='top', color=COL[k] if k != 'mean-effect' else '#636363', fontweight='bold')
letter(c, 'F', -0.42, 1.04)
save(fig, 'Fig_S_mass')
print('saved Fig_S_classification and Fig_S_mass')
