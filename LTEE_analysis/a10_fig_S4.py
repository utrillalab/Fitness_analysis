"""Revised Supplementary Figure 4: loss-of-function (LoF) mutations in the 12 LTEE populations vs glucose fitness.

Inputs: outputs of `a09_ltee.py plus|minus|all` (run those first).
  A  glucose fitness of all genes with coding mutations (any type, 12 populations) vs all genes
  B  % of genes per fitness category: genome background vs LoF genes
  C  robustness: % mean-effect (Wilson 95% CI) in each subset, with the genome background
  D  transporters with LoF mutations, by category
  E  transcriptional regulators with LoF mutations, by category
Category colours were checked for colour-vision deficiency (OKLab, Machado 2009 simulation): min pairwise
Delta E 13.7 (protan/deutan/tritan), 15.7 normal vision. Every bar is labelled, so colour is never the only cue.
Output: Output/Fig_S4.{svg,pdf,png}
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from common import *

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7.5, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
COL = {'essential': '#b2182b', 'important': '#e08214', 'mean-effect': '#969696', 'enhancing': '#2166ac'}
LAB = {'essential': 'Essential\n(≤ −2)', 'important': 'Important\n(−2 to −0.5)', 'mean-effect': 'Mean-effect\n(−0.5 to 0.5)',
       'enhancing': 'Enhancing\n(> 0.5)'}
BG = '#dcdcdc'
INK, INK2 = '#222222', '#666666'

# ---- data -------------------------------------------------------------------------------------------
fit, names = load_fitness()
glc = fit['Glucose']
eco = pd.read_csv(os.path.join(REPO, 'Input_Data', 'All-genes-of-E.-coli-K-12-substr-MG1655.txt'), sep='\t')
eco = eco.dropna(subset=['Accession-1', 'Left-End-Position', 'Right-End-Position'])   # same universe as a09_ltee.py
eco = eco[eco['Accession-1'].astype(str).str.match(r'^b\d{4}$')].drop_duplicates('Accession-1')
universe = [b for b in eco['Accession-1'] if b in glc.index and pd.notna(glc[b])]
bg = pd.Series(classify(glc.reindex(universe)), index=universe)

genes = pd.read_csv(f'{OUT}/a09_all_ltee_lof_genes.csv')
genes = genes[genes.in_universe]
tests = {k: pd.read_csv(f'{OUT}/a09{s}_ltee_lof_tests.csv') for k, s in [('plus', ''), ('minus', '_minus'), ('all', '_all')]}


def row(k, setname):
    t = tests[k]
    return t[(t.set == setname) & (t.categories == 'fitness only')].iloc[0]


def wilson(x, n, z=1.96):
    p = x / n; d = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / d
    return 100 * (c - h), 100 * (c + h)


def fmt_p(p):
    if p <= 1e-4:
        return r'$p < 10^{-4}$'
    if p >= 0.001:
        return f'$p = {p:.2g}$'
    m, e = f'{p:.0e}'.split('e')
    return rf'$p = {m}\times10^{{{int(e)}}}$'


# ---- figure -----------------------------------------------------------------------------------------
from matplotlib.patches import Patch
fig = plt.figure(figsize=(7.2, 8.2))
# three rows: A | B ; C (full width) ; D | E   (positions in figure fractions)
axA = fig.add_axes([0.085, 0.705, 0.385, 0.255])
axB = fig.add_axes([0.590, 0.705, 0.395, 0.255])
axC = fig.add_axes([0.380, 0.365, 0.400, 0.245])
axD = fig.add_axes([0.085, 0.065, 0.385, 0.200])
axE = fig.add_axes([0.590, 0.065, 0.395, 0.200])


def letter(ax, s, x=None):
    # panel letters aligned on two columns in figure coordinates
    fx = 0.012 if x is None else x
    y = ax.get_position().y1 + 0.012
    fig.text(fx, y, s, fontsize=10, fontweight='bold', va='bottom', ha='left', color=INK)


# A: all genes with coding mutations (any type) in the 12 populations, parsed with a09_ltee.py (sections 1-2)
import sys, io, contextlib
from scipy import stats
_argv = sys.argv; sys.argv = ['a09_ltee.py', 'all']
_ns = {}
with contextlib.redirect_stdout(io.StringIO()):
    exec(open(os.path.join(HERE, 'a09_ltee.py'), encoding='utf-8').read().split('# 3. hypermutator check')[0], _ns)
sys.argv = _argv
_coding = _ns['d'][_ns['d']['class'].isin(_ns['LOF'] | {'missense', 'synonymous', 'inframe_indel'})]
mutated = sorted({b for bs in _coding.bnums for b in bs} & set(universe))
me_mut = (bg.reindex(mutated) == 'mean-effect').mean() * 100
p_mut = stats.hypergeom.sf(int((bg.reindex(mutated) == 'mean-effect').sum()) - 1, len(universe),
                           int((bg == 'mean-effect').sum()), len(mutated))
print(f'Panel A: {len(mutated)} genes with coding mutations; mean-effect {me_mut:.1f}% vs {100 * (bg == "mean-effect").mean():.1f}% '
      f'of all genes (hypergeometric p = {p_mut:.2f})')
pd.DataFrame({'bnumber': mutated, 'glucose_fitness': glc.reindex(mutated).values,
              'category': bg.reindex(mutated).values}).to_csv(f'{OUT}/a10_all_mutated_genes_12pops.csv', index=False)
bins = np.arange(-4, 1.6, 0.1)
fa = glc.reindex(universe).clip(-4, 1.5); fm = glc.reindex(mutated).clip(-4, 1.5)
axA.hist(fa, bins, density=True, color=BG, label=f'All genes (n = {len(universe):,})')
axA.hist(fm, bins, density=True, histtype='step', color=INK, lw=1.0, label=f'Genes with coding mutations (n = {len(mutated)})')
for v in (-2, -0.5, 0.5):
    axA.axvline(v, color=INK2, lw=0.5, ls=(0, (2, 2)), zorder=0)
axA.text(0.02, 0.62, f'Mean-effect: {me_mut:.1f}% of mutated genes\nvs {100 * (bg == "mean-effect").mean():.1f}% of all genes ($p = {p_mut:.2f}$)',
         transform=axA.transAxes, ha='left', va='top', fontsize=6.3, color=INK,
         bbox=dict(facecolor='white', edgecolor='none', pad=1.5), zorder=5)
axA.set_title('All mutated genes', fontsize=7, loc='left', color=INK, pad=4)
axA.set_xlim(-4, 1.5); axA.set_xlabel('Fitness in glucose (values below −4 shown at −4)'); axA.set_ylabel('Density')
leg = axA.legend(frameon=True, fontsize=6.5, loc='upper left', handlelength=1.4, borderpad=0.3, framealpha=1)
leg.get_frame().set_linewidth(0); leg.set_zorder(5)
letter(axA, 'A')

# B: categories, background vs LoF (12 populations)
r_all = row('all', 'all LoF')
w = 0.38
for i, k in enumerate(CATS):
    pb = 100 * (bg == k).mean(); pl = r_all[f'{k}_%']; nl = int(r_all[f'{k}_n'])
    axB.bar(i - w / 2 - 0.01, pb, w, color=BG, linewidth=0)
    axB.bar(i + w / 2 + 0.01, pl, w, color=COL[k], linewidth=0)
    axB.text(i - w / 2 - 0.01, pb + 1.5, f'{pb:.1f}%', ha='center', va='bottom', fontsize=6, color=INK2)
    axB.text(i + w / 2 + 0.01, pl + 1.5, f'{pl:.1f}%\n(n = {nl})', ha='center', va='bottom', fontsize=6, color=INK)
axB.set_xticks(range(4), [LAB[k] for k in CATS], fontsize=6.5); axB.set_ylim(0, 115); axB.set_yticks([0, 25, 50, 75, 100])
axB.set_ylabel('% of genes')
axB.legend(handles=[Patch(color=BG, label=f'All genes (n = {len(universe):,})'),
                    Patch(facecolor='white', edgecolor=INK, lw=0.6, label=f'Genes with LoF (n = {int(r_all.genes_in_universe)});\ncoloured by category')],
           frameon=False, fontsize=6.5, loc='upper left', handlelength=1.2)
exp_ess = r_all['essential_background_%'] * r_all.genes_in_universe / 100
axB.text(0.02, 0.62, f'Mean-effect enriched: {fmt_p(r_all["mean-effect_p_perm_enriched"])}\n'
                     f'Essential: {int(r_all.essential_n)} observed, {exp_ess:.1f} expected',
         transform=axB.transAxes, ha='left', va='top', fontsize=6.3, color=INK)
letter(axB, 'B', 0.505)

# C: robustness
rows = [('Ara+ populations', row('plus', 'all LoF'), True),
        ('Ara− populations', row('minus', 'all LoF'), True),
        ('All 12 populations', r_all, True),
        ('12 populations, without IS-mutator clones', row('all', 'LoF without IS-mutator clones'), True),
        ('12 populations, point mutations only', row('all', 'point-like LoF (no multi-gene deletions)'), True),
        ('12 populations, LoF in ≥ 2 populations', row('all', 'LoF in >=2 populations (parallel)'), True),
        ('Comparison: genes with missense mutations', row('all', 'comparator: missense'), False),
        ('Comparison: genes with synonymous mutations', row('all', 'comparator: synonymous'), False)]
base = r_all['mean-effect_background_%']
axC.axvline(base, color=INK2, lw=0.6, ls=(0, (2, 2)))
axC.text(base, len(rows) - 0.45, f'All genes, {base:.1f}%', ha='center', va='bottom', fontsize=6, color=INK2,
          bbox=dict(facecolor='white', edgecolor='none', pad=1))
for j, (lab, r, is_lof) in enumerate(rows):
    y = len(rows) - 1 - j
    n = int(r.genes_in_universe); xm = int(r['mean-effect_n']); lo, hi = wilson(xm, n)
    axC.plot([lo, hi], [y, y], color=INK, lw=0.8, solid_capstyle='butt', zorder=2)
    axC.plot(100 * xm / n, y, 'o', ms=4.5, mfc=COL['mean-effect'] if is_lof else 'white', mec=INK, mew=0.8, zorder=3)
    pp = r['mean-effect_p_perm_enriched'] if pd.notna(r.get('mean-effect_p_perm_enriched', np.nan)) else r['mean-effect_p_hyper_enriched']
    axC.text(1.03, y, f'n = {n}' + (f',  {fmt_p(pp)}' if is_lof else ''), transform=axC.get_yaxis_transform(),
             va='center', fontsize=6.3, color=INK)
axC.set_yticks(range(len(rows)), [r[0] for r in rows][::-1], fontsize=6.5)
axC.set_xlim(70, 100); axC.set_xticks([70, 75, 80, 85, 90, 95, 100]); axC.set_ylim(-0.6, len(rows) + 0.3)
axC.set_xlabel('% mean-effect genes (95% CI)')
axC.spines['left'].set_visible(False); axC.tick_params(axis='y', length=0)
letter(axC, 'C')

# D/E: transporters and regulators
for ax_, col, title, let, lx in [(axD, 'transporter (EcoCyc product)', 'Transporters', 'D', None),
                                 (axE, 'RegulonDB_regulator', 'Transcriptional regulators', 'E', 0.505)]:
    s = genes[genes[col]]
    cnt = [int((s.category == k).sum()) for k in CATS]
    ax_.bar(range(4), cnt, 0.6, color=[COL[k] for k in CATS], linewidth=0)
    top = max(cnt) * 1.25
    for i, (k, v) in enumerate(zip(CATS, cnt)):
        names_k = s[s.category == k].gene.tolist()
        ax_.text(i, v + top * 0.02, str(v), ha='center', va='bottom', fontsize=6.3, color=INK)
        if 0 < v <= 3:
            ax_.text(i, v + top * 0.09, ', '.join(names_k), ha='center', va='bottom', fontsize=6.3, color=INK,
                     fontstyle='italic')
    ax_.set_xticks(range(4), [LAB[k] for k in CATS], fontsize=6.3)
    ax_.set_ylim(0, top); ax_.set_ylabel('Genes with LoF')
    ax_.set_title(f'{title} (n = {len(s)})', fontsize=7, color=INK, pad=4, loc='left')
    letter(ax_, let, lx)

for ext in ('svg', 'pdf'):
    fig.savefig(f'{OUT}/Fig_S4.{ext}')
fig.savefig(f'{OUT}/Fig_S4.png', dpi=300)
print('saved', f'{OUT}/Fig_S4.*', '| universe', len(universe), '| LoF genes', len(genes))
