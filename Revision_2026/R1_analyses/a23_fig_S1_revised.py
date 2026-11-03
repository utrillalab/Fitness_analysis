"""Revised Supplementary Figure S1 (draft layout for review).

New panels drawn here:
  * protein mass by fitness category in each condition (replaces the original S1 box plots, which used the misaligned
    Table S1 fitness): box plots of log10 protein mass, with Cliff's delta vs mean-effect and Holm-adjusted
    Kruskal-Wallis significance (R1-M4b; output/a02_M4b_effect_sizes.csv);
  * condition specificity of each category (R1-M8; output/a03_M8_condition_specificity.csv): mean pairwise Jaccard
    index between the 11 conditions vs a size-matched random null, and number of conditions in which each gene keeps
    its category.
The draft layout also places panels made by other scripts (as images): fitness distributions (a01, R1-M2), t statistic
(a05, R2-M2a), proteins without fitness data (a13, R1-M1), mass quintiles (a18, R1-M4c).
Output: output/Fig_S1_condition_specificity.{svg,pdf,png}, output/Fig_S1_mass_boxplots.{svg,pdf,png},
        output/Fig_S1_revised_draft_p1.png, output/Fig_S1_revised_draft_p2.png
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from common import *

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
COL = {'essential': '#b2182b', 'important': '#e08214', 'mean-effect': '#969696', 'enhancing': '#2166ac'}
LAB = {'essential': 'Essential', 'important': 'Important', 'mean-effect': 'Mean-effect', 'enhancing': 'Enhancing'}
ORDER = ['mean-effect', 'important', 'essential', 'enhancing']
INK, INK2 = '#222222', '#666666'


def save(fig, name):
    for ext in ('svg', 'pdf'):
        fig.savefig(f'{OUT}/{name}.{ext}')
    fig.savefig(f'{OUT}/{name}.png', dpi=300)


def letter(ax, s, x=-0.18, y=1.08):
    ax.text(x, y, s.lower(), transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom')  # npj: lower-case bold


# ---- condition specificity (R1-M8) ----------------------------------------------------------------
fit, _ = load_fitness()
cat11 = pd.DataFrame({c: classify(fit[c]) for c in CONDS}, index=fit.index)
spec = pd.read_csv(f'{OUT}/a03_M8_condition_specificity.csv').set_index('category')


def panel_jaccard(ax):
    ks = ['mean-effect', 'essential', 'important', 'enhancing']
    x = np.arange(len(ks))
    ax.bar(x, spec.loc[ks, 'mean_pairwise_Jaccard_between_conditions'], 0.6, color=[COL[k] for k in ks])
    ax.scatter(x, spec.loc[ks, 'Jaccard_null_mean'], marker='_', s=180, color='black', lw=1.2, zorder=3,
               label='Random sets of\nthe same size')
    for i, k in enumerate(ks):
        ax.text(i, spec.loc[k, 'mean_pairwise_Jaccard_between_conditions'] + 0.02,
                f"{spec.loc[k, 'mean_pairwise_Jaccard_between_conditions']:.2f}", ha='center', va='bottom', fontsize=7)
    ax.set_xticks(x, [LAB[k] for k in ks], rotation=30, ha='right')
    ax.set_ylim(0, 1); ax.set_ylabel('Gene-set overlap between conditions\n(mean pairwise Jaccard index)')
    ax.legend(frameon=False, fontsize=6.5, loc='upper right', handlelength=1)


def panel_nconds(ax):
    n = np.arange(1, 12)
    for k in ['essential', 'important', 'enhancing', 'mean-effect']:
        cnt = (cat11 == k).sum(axis=1); cnt = cnt[cnt > 0]
        frac = np.array([(cnt == i).mean() * 100 for i in n])
        ax.plot(n, frac, '-o', color=COL[k], ms=2.5, lw=1.1,
                label=f'{LAB[k]} (median {int(spec.loc[k, "median_conditions_per_gene (of 11)"])})')
    ax.set_xticks(n); ax.set_xlim(0.5, 11.5)
    ax.set_xlabel('Conditions in which the gene is in the category'); ax.set_ylabel('% of genes in the category')
    ax.legend(frameon=False, fontsize=6.5, loc='upper center', handlelength=1.2)


fig, (a, b) = plt.subplots(1, 2, figsize=(5.6, 2.4), gridspec_kw={'width_ratios': [0.8, 1.2], 'wspace': 0.45})
panel_jaccard(a); panel_nconds(b); letter(a, 'X', -0.42); letter(b, 'Y', -0.22)
fig.subplots_adjust(left=0.13, right=0.98, bottom=0.25, top=0.86)
save(fig, 'Fig_S1_condition_specificity')

# ---- mass by category per condition (replaces the original S1 box plots) ---------------------------
df = long_table(); df = df[df.mass_fg > 0].assign(logm=lambda x: np.log10(x.mass_fg))
eff = pd.read_csv(f'{OUT}/a02_M4b_effect_sizes.csv').set_index(['condition', 'category'])
CORDER = ['Glucose', 'Acetate', 'Glucosamine', 'Glycerol', 'Pyruvate', 'Xylose', 'Mannose', 'Galactose', 'Fructose',
          'Succinate', 'LB']


def stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'n.s.'


def boxes(axs):
    for ax, c in zip(axs, CORDER):
        d = df[df.condition == c]
        data = [d.loc[d.category == k, 'logm'].values for k in ORDER]
        bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color='black', lw=0.8),
                        whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6), boxprops=dict(lw=0.6))
        for patch, k in zip(bp['boxes'], ORDER):
            patch.set_facecolor(COL[k]); patch.set_alpha(0.85)
        for i, k in enumerate(ORDER[1:], start=2):
            if (c, k) in eff.index and len(data[i - 1]) >= 3:
                r = eff.loc[(c, k)]
                ax.text(i, 1.6, f"{stars(r.p_holm)}\n{r.cliffs_delta:.2f}", ha='center', va='bottom', fontsize=5.3,
                        color=INK)
        ax.set_ylim(-6.5, 2.6); ax.set_title(c, fontsize=7.5, pad=12)
        ax.set_xticks(range(1, 5), ['ME', 'Imp', 'Ess', 'Enh'], fontsize=6.5)
        ax.text(1, -6.3, f'n={len(data[0])}', ha='center', fontsize=5, color=INK2)


fig, axs = plt.subplots(2, 6, figsize=(7.2, 3.6), sharey=True)
axs = axs.ravel()
boxes(axs[:11])
axs[11].axis('off')
axs[11].text(-0.1, 0.5, 'ME, mean-effect\nImp, important\nEss, essential\nEnh, enhancing\n\nNumbers: Cliff\'s delta\n'
             'vs mean-effect\n* Holm-adjusted p < 0.05\n** < 0.01, *** < 0.001\n(Kruskal-Wallis)', fontsize=6.3, va='center', color=INK)
for ax in axs[::6]:
    ax.set_ylabel('log$_{10}$ protein mass (fg)')
fig.subplots_adjust(left=0.08, right=0.99, bottom=0.08, top=0.9, hspace=0.55, wspace=0.12)
save(fig, 'Fig_S1_mass_boxplots')


# ---- draft layout (2 pages, existing panels placed as images) ----------------------------------------
def place(fig, rect, path, lab, title=None):
    ax = fig.add_axes(rect); ax.imshow(mpimg.imread(f'{OUT}/{path}')); ax.axis('off')
    fig.text(rect[0], rect[1] + rect[3] + 0.003, lab, fontsize=12, fontweight='bold', va='bottom')
    if title:
        fig.text(rect[0] + 0.03, rect[1] + rect[3] + 0.005, title, fontsize=7, color=INK2, va='bottom')


p1 = plt.figure(figsize=(8.27, 11.69))
p1.text(0.5, 0.985, 'Supplementary Figure S1 (draft layout, page 1): fitness distributions and classification',
        ha='center', fontsize=9, fontweight='bold')
place(p1, [0.04, 0.60, 0.92, 0.36], 'a01_fitness_distributions.png', 'A', 'Fitness distribution per condition (R1-M2)')
place(p1, [0.04, 0.33, 0.92, 0.23], 'a05_t_statistic.png', 'B', 't-like statistic (R2-M2a)')
place(p1, [0.10, 0.02, 0.80, 0.28], 'a13_missing_fitness.png', 'C', 'Proteins without fitness data (R1-M1)')
p1.savefig(f'{OUT}/Fig_S1_revised_draft_p1.png', dpi=200)

p2 = plt.figure(figsize=(8.27, 11.69))
p2.text(0.5, 0.985, 'Supplementary Figure S1 (draft layout, page 2): protein mass and fitness categories',
        ha='center', fontsize=9, fontweight='bold')
place(p2, [0.03, 0.56, 0.94, 0.40], 'Fig_S1_mass_boxplots.png', 'D', 'Protein mass by category, per condition (original S1, corrected; R1-M4b)')
place(p2, [0.03, 0.30, 0.94, 0.23], 'Fig_S1_mass_quintiles.png', 'E', 'Mass distribution and quintile enrichment (R1-M4c)')
place(p2, [0.12, 0.02, 0.76, 0.25], 'Fig_S1_condition_specificity.png', 'F', 'Condition specificity (R1-M8)')
p2.savefig(f'{OUT}/Fig_S1_revised_draft_p2.png', dpi=200)
print('saved Fig_S1_condition_specificity, Fig_S1_mass_boxplots, Fig_S1_revised_draft_p1/p2')
