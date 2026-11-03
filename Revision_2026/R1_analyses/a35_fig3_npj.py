"""Fig. 3 (GO dot plot) redrawn in Python in the npj style (28 Oct 2026; R is not installed on this computer).

Same data as a15c_fig3_dotplot.R, read from its output a15c_fig3_terms.csv (run "corrected"): one dot per GO term and
category; colour = value (the -log10 p-value of the term, mean over the carbon sources, min-max scaled to 0-10 in each
category, as in the R script), size = freq (fraction of the genes of the term, also min-max scaled). Terms in reverse
alphabetical order from top to bottom and dot area scaled as in ggplot2 (scale_size; range 1-4.5 mm instead of 1-6 mm, for the
smaller print size), viridis colours.
Changes: Arial 7 pt at print size (150 mm wide), colour legend "Scaled -log10(p value)" (0-10 within each
category, as stated in the figure legend) with a proper minus sign and italic p.
Output: output/Fig03_dotplot_corrected.{pdf,svg,png}
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from npj_style import *
OUT = os.path.join(HERE, 'output')

d = pd.read_csv(os.path.join(OUT, 'a15c_fig3_terms.csv'))
d = d[d.run == 'corrected'].copy()
d['variable'] = d.variable.replace({'No effect': 'Mean-effect'})
CATS = ['Essential', 'Important', 'Mean-effect']
terms = sorted(d.term.unique(), key=str.lower)          # ggplot: first level at the bottom
y = {t: i for i, t in enumerate(terms)}
x = {c: i for i, c in enumerate(CATS)}

PT = 72 / 25.4                                            # points per mm
lo, hi = d.freq.min(), d.freq.max()


def area(f):                                              # ggplot2 scale_size: area linear in the value, 1-4.5 mm diameter
    r = np.sqrt((np.asarray(f) - lo) / (hi - lo) * (4.5 ** 2 - 1 ** 2) + 1 ** 2)
    return (r * PT) ** 2


fig = plt.figure(figsize=(150 * MM, 120 * MM))
ax = fig.add_axes([0.44, 0.08, 0.34, 0.90])
for t in terms:
    ax.axhline(y[t], color='#ebebeb', lw=0.5, zorder=0)
for c in CATS:
    ax.axvline(x[c], color='#ebebeb', lw=0.5, zorder=0)
sc = ax.scatter(d.variable.map(x), d.term.map(y), s=area(d.freq), c=d.value, cmap='viridis', lw=0, zorder=2)
ax.set_yticks(range(len(terms)), terms); ax.set_xticks(range(3), CATS)
ax.set_xlim(-0.6, 2.6); ax.set_ylim(-0.7, len(terms) - 0.3)
ax.set_xlabel('Categories')
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)

cax = fig.add_axes([0.82, 0.60, 0.02, 0.16])
cb = fig.colorbar(sc, cax=cax); cb.outline.set_visible(False); cb.ax.tick_params(length=2); cb.set_ticks([4, 6, 8, 10])
cax.set_title('Scaled\n−log$_{10}$($\\it{p}$ value)', loc='left', pad=6)
br = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
h = [plt.scatter([], [], s=area(b), color='black', lw=0) for b in br]
fig.legend(h, [f'{b:.1f}' for b in br], title='Fraction of genes\nin the GO term', loc='upper left',
           bbox_to_anchor=(0.81, 0.50), labelspacing=0.6, alignment='left', borderaxespad=0)
save(fig, os.path.join(OUT, 'Fig03_dotplot_corrected'))
print(len(d), 'dots,', len(terms), 'terms')
