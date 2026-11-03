"""New Supplementary Figure S1 panels (R1-M4c): protein mass by fitness category and category representation per
protein-mass quintile.

  A  distribution of log10 protein mass per category, 11 conditions pooled (kernel density); shaded bands = median
     quintile boundaries across conditions (the test in B uses each condition's own quintiles).
  B  observed/expected genes per category and quintile, pooled over the 11 conditions, and the number of conditions
     with significant over- (↑) or under-representation (↓) (hypergeometric, BH over all tests, p < 0.05).
Inputs: common.long_table() and output/a02_M4c_pooled.csv (a02_mass_by_category.py). Proteins with zero mass are
excluded, as in a02. Style and category colours as Supplementary Figure 4 (colour-vision-deficiency checked).
Output: output/Fig_S1_mass_quintiles.{svg,pdf,png} (A + B) and output/Fig_S1_mass_quintiles_B.{svg,pdf,png} (B only)
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import gaussian_kde
from common import *

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
COL = {'essential': '#b2182b', 'important': '#e08214', 'mean-effect': '#636363', 'enhancing': '#2166ac'}
LAB = {'essential': 'Essential', 'important': 'Important', 'mean-effect': 'Mean-effect', 'enhancing': 'Enhancing'}
BINS = ['Q1 (lightest)', 'Q2', 'Q3', 'Q4', 'Q5 (heaviest)']
NC = 11
XMIN, XMAX = -7, 1.5
INK = '#222222'

df = long_table()
df = df[df.mass_fg > 0].copy()
df['logm'] = np.log10(df.mass_fg)
edges = df.groupby('condition').logm.quantile([0.2, 0.4, 0.6, 0.8]).unstack().median().values
p = pd.read_csv(f'{OUT}/a02_M4c_pooled.csv').set_index(['category', 'mass_bin'])
R = np.array([[p.loc[(k, b), 'obs/exp'] for b in BINS] for k in CATS])
L = np.log2(R)


def letter(ax, s, x=-0.16):
    ax.text(x, 1.06, s, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')


def panel_A(ax):
    bounds = [XMIN, *edges, XMAX]
    for j in range(5):
        if j % 2 == 0:
            ax.axvspan(bounds[j], bounds[j + 1], color='#f0f0f0', lw=0, zorder=0)
        ax.text((bounds[j] + bounds[j + 1]) / 2, 1.0, f'Q{j + 1}',
                transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=6.5, color='#666666')
    xs = np.linspace(XMIN, XMAX, 400)
    for k in CATS:
        v = df.loc[df.category == k, 'logm']
        ax.plot(xs, gaussian_kde(v, bw_method=0.15)(xs), color=COL[k], lw=1.3,
                label=f'{LAB[k]} (n = {len(v):,})', zorder=3)
    ax.set_xlim(XMIN, XMAX); ax.set_ylim(0, None)
    ax.set_xlabel('log$_{10}$ protein mass (fg per cell)'); ax.set_ylabel('Density')
    ax.legend(frameon=False, fontsize=6.5, loc='upper left', handlelength=1.5, borderaxespad=0.2)


def panel_B(ax, fig):
    im = ax.imshow(L, cmap='RdBu_r', norm=TwoSlopeNorm(0, -2, 2), aspect='auto')
    for i, k in enumerate(CATS):
        for j, b in enumerate(BINS):
            up, dn = int(p.loc[(k, b), 'n_cond_over']), int(p.loc[(k, b), 'n_cond_under'])
            sig = f'↑ {up}/{NC}' if up >= dn and up else (f'↓ {dn}/{NC}' if dn else 'n.s.')
            ink = 'white' if abs(L[i, j]) > 1.1 else INK
            ax.text(j, i - 0.14, f'{R[i, j]:.2f}×', ha='center', va='center', fontsize=7, color=ink, fontweight='bold')
            ax.text(j, i + 0.2, sig, ha='center', va='center', fontsize=6, color=ink)
    ax.set_xticks(range(5), ['Q1\n(lightest)', 'Q2', 'Q3', 'Q4', 'Q5\n(heaviest)'])
    ax.set_yticks(range(4), [LAB[k] for k in CATS])
    for t, k in zip(ax.get_yticklabels(), CATS):
        t.set_color(COL[k]); t.set_fontweight('bold')
    ax.tick_params(length=0)
    ax.set_xlabel('Protein-mass quintile (per condition)')
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks(np.arange(-0.5, 5), minor=True); ax.set_yticks(np.arange(-0.5, 4), minor=True)
    ax.grid(which='minor', color='white', lw=1.2); ax.tick_params(which='minor', length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03, ticks=[-2, -1, 0, 1, 2])
    cb.ax.set_yticklabels(['0.25×', '0.5×', '1×', '2×', '4×'], fontsize=6.5)
    cb.set_label('Observed / expected', fontsize=7); cb.outline.set_linewidth(0.4); cb.ax.tick_params(length=2, width=0.4)


fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.6), gridspec_kw={'width_ratios': [1, 1.25], 'wspace': 0.45})
panel_A(a); panel_B(b, fig)
letter(a, 'A', -0.2); letter(b, 'B', -0.3)
fig.subplots_adjust(left=0.08, right=0.92, bottom=0.2, top=0.88)
for ext in ('svg', 'pdf'):
    fig.savefig(f'{OUT}/Fig_S1_mass_quintiles.{ext}')
fig.savefig(f'{OUT}/Fig_S1_mass_quintiles.png', dpi=300)

fig2, ax2 = plt.subplots(figsize=(4.2, 2.3))
panel_B(ax2, fig2); fig2.tight_layout()
for ext in ('svg', 'pdf'):
    fig2.savefig(f'{OUT}/Fig_S1_mass_quintiles_B.{ext}')
fig2.savefig(f'{OUT}/Fig_S1_mass_quintiles_B.png', dpi=300)
print('quintile edges (median across conditions, log10 fg):', edges.round(2),
      '| proteins below 10^-7 fg not shown in A:', int((df.logm < XMIN).sum()), 'of', len(df))
