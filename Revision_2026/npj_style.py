"""Common style of the main figures for npj Systems Biology and Applications (28 Oct 2026).

npj artwork guide: Arial/Helvetica, the same font in every figure, ~8 pt at print size (7 pt here, 5-8 pt allowed),
bold lower-case panel letters, white background, colour-blind-safe colours, vector files.
Every figure is drawn at its print width (DOUBLE = 180 mm, two columns; SINGLE = 88 mm) so font sizes are the printed
ones. Category colours are those of Fig. 4 (essential red, important orange, mean-effect grey, enhancing blue).
Usage: from npj_style import *  (root scripts) or sys.path.insert(0, <Revision_2026>) first (R1_analyses scripts).
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MM = 1 / 25.4
DOUBLE, SINGLE = 180 * MM, 88 * MM
FS = 7          # body text, pt at print size
FS_LETTER = 8   # panel letters, bold
PAL = {'essential': '#b2182b', 'important': '#e08214', 'mean-effect': '#969696', 'enhancing': '#2166ac'}
LABEL = {'essential': 'Essential', 'important': 'Important', 'mean-effect': 'Mean-effect', 'enhancing': 'Enhancing'}
GREYS = ('#d9d9d9', '#969696', '#404040')   # light / mid / dark grey (subsets of the mean-effect genes, Fig. 1d)
BAR = '#5f9ea0'                             # single-series bars (Figs. 5 and 6)

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': FS, 'axes.titlesize': FS, 'axes.labelsize': FS, 'xtick.labelsize': FS,
    'ytick.labelsize': FS, 'legend.fontsize': FS, 'figure.titlesize': FS, 'axes.titleweight': 'normal',
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.minor.size': 1.5,
    'ytick.minor.size': 1.5, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.unicode_minus': True,
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold'})


def letter(fig, ax, s, dx=-0.02, dy=0.0):
    """Bold lower-case panel letter at the top-left corner of the axes' bounding box (incl. tick labels)."""
    fig.canvas.draw()
    bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bb.x0 + dx, bb.y1 + dy, s, fontsize=FS_LETTER, fontweight='bold', ha='left', va='bottom')


def save(fig, path_noext, exts=('pdf', 'svg', 'png')):
    """Vector files plus a 600-dpi RGB PNG (the PNG goes into the Word manuscript)."""
    from PIL import Image
    for e in exts:
        fig.savefig(f'{path_noext}.{e}', dpi=600)
    if 'png' in exts:
        Image.open(f'{path_noext}.png').convert('RGB').save(f'{path_noext}.png', dpi=(600, 600))
