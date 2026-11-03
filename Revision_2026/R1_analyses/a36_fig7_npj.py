"""Fig. 7 (J. A. Farias, SVG of 22 Oct 2026, external/Figure7_JAF.svg) adapted to the npj style (28 Oct 2026,
approved by the user). The SVG is matplotlib output with the text drawn as Helvetica Neue glyph outlines; each text
block carries its string in a comment and its position and size in its transform, so every block is replaced by real
text in Arial at the same position and size (textLength keeps the original width, so the layout does not move).
  - panel letters a-e in bold, 8 pt at the printed width (180 mm);
  - gene names in italics (ydcI, yceQ, ydhL, yobF, yjjY, tsaMBCD); "SEQUENCE"/"STRUCTURE" -> "Sequence"/"Structure";
  - confidence colours (panels a, b) colour-blind safe: high confidence dark green -> blue, lower confidence dark yellow
    -> amber, no hits wine red -> grey (the domain and DNA colours of c-e are kept: they match the rendered structures);
  - no frames around panels d and e.
Outputs (R1_analyses/output/): Fig07_npj.svg, Fig07_npj.png (600 dpi at 180 mm, rendered with Edge, RGB)
Run from Revision_2026: python R1_analyses/a36_fig7_npj.py
"""
import os, re, html, subprocess, tempfile
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(os.path.dirname(HERE), 'external', 'Figure7_JAF.svg')
OUT = os.path.join(HERE, 'output')
s = open(SRC, encoding='utf8').read()
W_PT = float(re.search(r'viewBox="0 0 ([\d.]+) ', s).group(1))
K = 180 / 25.4 * 72 / W_PT                     # print scale (180 mm wide)
GENES = ['ydcI', 'yceQ', 'ydhL', 'yobF', 'yjjY', 'tsaMBCD']
RENAME = {'SEQUENCE': 'Sequence', 'STRUCTURE': 'Structure'}
COLOURS = {'#1e6b47': '#1f5fa8', '#9a7a0a': '#d08c00', '#73243a': '#7f7f7f'}

# glyph widths (max x of each outline; outlines are stored at scale(1/64)) to measure the original width of each line
GW = {g: max([float(v) for v in re.findall(r'(-?[\d.]+) -?[\d.]+', d)] or [0]) / 64
      for g, d in re.findall(r'<path id="([^"]+)" d="([^"]*)"', s)}


def spans(t):
    t = html.escape(t, quote=False)
    for g in GENES:
        t = re.sub(rf'\b{g}\b', f'<tspan font-style="italic">{g}</tspan>', t)
    return t


def repl(m):
    txt, style, tr, body = m.groups()
    sc = float(re.search(r'scale\(([\d.]+)', tr).group(1))
    fs = sc * 100                               # matplotlib text-as-path: glyphs at 100 units per point of font size
    uses = re.findall(r'<use xlink:href="#([^"]+)"(?: x="([-\d.]+)"| transform="translate\(([-\d.]+) )?', body)
    width = max((float(u[1] or u[2] or 0) + GW.get(u[0], 0)) * sc for u in uses) if uses else 0
    fill = re.search(r'fill: (#[0-9a-fA-F]{6})', style or '')
    fill = fill.group(1) if fill else '#000000'
    pos = re.match(r'translate\(([-\d.]+) ([-\d.]+)\)( rotate\([-\d.]+\))?', tr)
    x, y, rot = float(pos.group(1)), float(pos.group(2)), pos.group(3) or ''
    attrs = f'font-family="Arial, Helvetica, sans-serif" font-size="{fs:.3f}" fill="{fill}"'
    if re.fullmatch(r'[a-e]', txt) and fs > 15:           # panel letters: bold, 8 pt at print size
        return (f'<!-- {txt} --><text x="{x:.3f}" y="{y:.3f}" font-family="Arial, Helvetica, sans-serif" '
                f'font-size="{8 / K:.3f}" fill="{fill}" font-weight="bold">{txt}</text>')
    if txt in RENAME:                                      # new string: keep the centre, natural width
        return (f'<!-- {RENAME[txt]} --><text x="{x + width / 2:.3f}" y="{y:.3f}" text-anchor="middle" '
                f'{attrs}>{RENAME[txt]}</text>')
    tl = f' textLength="{width:.3f}" lengthAdjust="spacingAndGlyphs"' if width > 0 and len(txt) > 1 else ''
    if rot:
        return (f'<!-- {txt} --><text transform="translate({x:.3f} {y:.3f}){rot}" x="0" y="0" {attrs}{tl}>'
                f'{spans(txt)}</text>')
    return f'<!-- {txt} --><text x="{x:.3f}" y="{y:.3f}" {attrs}{tl}>{spans(txt)}</text>'


# one match per text line (multi-line blocks hold several comment + glyph-group pairs)
LINE = re.compile(r'<!-- (.*?) -->\s*<g (?:style="([^"]*)" )?transform="([^"]*)">((?:(?!</g>).)*)</g>', re.S)
n0 = s.count('<!--')
s, n = LINE.subn(repl, s)
assert n == n0 == 80, (n, n0)
assert '#HelveticaNeue' not in s, 'glyph outlines left'   # (scatter markers are <use> too)
for a, b in COLOURS.items():
    s = s.replace(a, b)
for p in ('patch_2', 'patch_3'):                           # frames of panels d and e
    s = re.sub(rf'(<g id="{p}">\s*<path [^>]*style="fill: #ffffff); stroke: #b9bfc4; stroke-width: 0.8; stroke-linejoin: miter"',
               r'\1"', s, count=1)
assert s.count('stroke: #b9bfc4; stroke-width: 0.8; stroke-linejoin: miter') < 10
svg = os.path.join(OUT, 'Fig07_npj.svg')
open(svg, 'w', encoding='utf8').write(s)

# render with Edge (headless) at 600 dpi for a 180-mm print width
H_PT = float(re.search(r'viewBox="0 0 [\d.]+ ([\d.]+)"', s).group(1))
wpx, hpx = round(W_PT * 4 / 3), round(H_PT * 4 / 3)
f = 180 / 25.4 * 600 / wpx
page = os.path.join(tempfile.gettempdir(), 'fig7_render.html')
open(page, 'w', encoding='utf8').write(f'<html><body style="margin:0;background:#fff"><img src="file:///{svg.replace(os.sep, "/")}" '
                                       f'width="{wpx}" height="{hpx}"></body></html>')
png = os.path.join(OUT, 'Fig07_npj.png')
EDGE = os.environ.get('EDGE', r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe')  # any Chromium browser
subprocess.run([EDGE, '--headless', '--disable-gpu', '--hide-scrollbars', f'--force-device-scale-factor={f:.4f}',
                f'--window-size={wpx},{hpx}', f'--screenshot={png}', 'file:///' + page.replace(os.sep, '/')],
               check=True, capture_output=True, timeout=120)
im = Image.open(png).convert('RGB')
im.save(png, dpi=(600, 600))
print('text blocks replaced:', n, '| print scale', round(K, 3), '| png', im.size)
