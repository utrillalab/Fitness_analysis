"""New supplementary figure (provisional S3; R1-M3 and R1-m3): ME-model predictions vs fitness categories and proteomics.

  A  % of the genes of each fitness category present in the ME-model that have zero predicted translation flux, in each
     of the 10 carbon sources (points) and their median (bar); mean-effect vs other categories: Fisher's exact test.
  B  glucose: predicted translation flux vs measured protein mass for genes with non-zero flux, coloured by category;
     Spearman rho in glucose and range over the 10 conditions.
Inputs: output/a03_M3_me_status_by_category.csv, output/a03_m3_flux_correlations.csv (a03_me_specificity_threshold.py),
common.load_flux() and common.long_table().
Output: output/Fig_S_me_model.{svg,pdf,png}
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from common import *

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
COL = {'essential': '#b2182b', 'important': '#e08214', 'mean-effect': '#969696', 'enhancing': '#2166ac'}
LAB = {'essential': 'Essential', 'important': 'Important', 'mean-effect': 'Mean-effect', 'enhancing': 'Enhancing'}
INK, INK2 = '#222222', '#666666'
ORDER = ['essential', 'important', 'enhancing', 'mean-effect']

st = pd.read_csv(f'{OUT}/a03_M3_me_status_by_category.csv')
z = st[st.me_status == 'zero flux'].pivot(index='condition', columns='category', values='%_of_category_genes_in_ME')
nme = st[st.me_status != 'absent from ME'].groupby(['condition', 'category']).n_genes.sum().unstack()
cor = pd.read_csv(f'{OUT}/a03_m3_flux_correlations.csv').set_index('condition')

fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.42})
rng = np.random.default_rng(0)
for i, k in enumerate(ORDER):
    v = z[k].dropna()
    a.barh(i, v.median(), 0.6, color=COL[k], alpha=0.9)
    a.scatter(v, i + rng.uniform(-0.18, 0.18, len(v)), s=7, color='white', edgecolor=INK, lw=0.5, zorder=3)
    a.text(101, i, f'{v.median():.0f}%\n(n ≈ {int(nme[k].median())})', va='center', fontsize=6.5, color=INK)
a.set_yticks(range(4), [LAB[k] for k in ORDER]); a.set_xlim(0, 100)
a.set_xlabel('Genes in the ME-model with zero\npredicted translation flux (%)')
a.text(0.02, -0.38, f'Mean-effect vs other categories: odds ratio {cor.OR_zero_flux_mean_effect_vs_other.min():.1f}–'
       f'{cor.OR_zero_flux_mean_effect_vs_other.max():.1f}, p < 10$^{{-10}}$ in 10 of 10 conditions (Fisher)',
       transform=a.transAxes, fontsize=6.3, color=INK2)
a.text(-0.34, 1.04, 'a', transform=a.transAxes, fontsize=10, fontweight='bold')

flux = load_flux()['Glucose']
d = long_table(); d = d[d.condition == 'Glucose'].assign(flux=lambda x: x.bnum.map(flux))
d = d[(d.flux > 0) & (d.mass_fg > 0)]
for k in ['mean-effect', 'important', 'essential', 'enhancing']:
    s = d[d.category == k]
    b.scatter(np.log10(s.flux), np.log10(s.mass_fg), s=6, color=COL[k], lw=0, alpha=0.85, label=f'{LAB[k]} (n = {len(s)})')
rho = stats.spearmanr(d.flux, d.mass_fg).statistic
b.set_xlabel('log$_{10}$ predicted translation flux (ME-model)'); b.set_ylabel('log$_{10}$ measured protein mass (fg)')
b.text(0.02, 0.97, f'Glucose: Spearman ρ = {rho:.2f}\n10 carbon sources: ρ = {cor["spearman_flux_vs_mass (non-zero)"].min():.2f}–'
       f'{cor["spearman_flux_vs_mass (non-zero)"].max():.2f}, p < 10$^{{-20}}$', transform=b.transAxes, va='top', fontsize=6.5)
b.legend(frameon=False, fontsize=6.3, loc='lower right', markerscale=1.8, handletextpad=0.1)
b.text(-0.2, 1.04, 'b', transform=b.transAxes, fontsize=10, fontweight='bold')
fig.subplots_adjust(left=0.12, right=0.98, bottom=0.3, top=0.93)
for ext in ('svg', 'pdf'):
    fig.savefig(f'{OUT}/Fig_S_me_model.{ext}', dpi=600)
fig.savefig(f'{OUT}/Fig_S_me_model.png', dpi=300)
print(z[ORDER].median().round(1).to_dict(), 'rho glucose', round(rho, 3), 'n', len(d))
