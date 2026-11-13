"""MUBLO group drawn in the style of Kaminski+ 2018 Fig. (n_X/n_CO bars), as a separate panel whose y axis matches
kaminski2018_abundances.png pixel-for-pixel (same 1258 px height): 10^-5 at y=1154.2, 255.67 px/decade, frame 8-12 / 1152-1156.
MUBLO values from the MUBLO analysis (abundance_comparison.ecsv, 2026-10-28): LTE at T_ex = 11.6 K, 0.5" aperture,
CO from 13CO 3-2 at 12C/13C = 25; gray bars 12C/13C = 4-89 (X/CO scales as 25/R).  SiO is a 3-sigma limit."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
out = sys.argv[1]
plt.rcParams.update({"font.family": "Times New Roman", "mathtext.fontset": "stix"})
W, H, dpi = 480, 1258, 100
PT = 72.0 / dpi                                   # points per pixel
yb, dec = 1154.2, 255.67                          # 10^-5 row and px/decade in the K18 image
ytop = 10.0
lo, hi = -5.0, -5.0 + (yb - ytop) / dec
xl, xr = 42.0, 470.0                              # spine centers (K18 group pitch 428.5 px)
fig = plt.figure(figsize=(W / dpi, H / dpi), dpi=dpi)
ax = fig.add_axes([xl / W, (H - yb) / H, (xr - xl) / W, (yb - ytop) / H])
ax.set_xlim(0, xr - xl); ax.set_ylim(lo, hi)
for s in ax.spines.values(): s.set_linewidth(5 * PT)
ax.set_xticks([(xr - xl) / 2]); ax.set_xticklabels([])
ax.tick_params(axis="x", which="both", direction="in", top=True, bottom=False, length=47 * PT, width=3 * PT)
maj = np.arange(-5, 0); mnr = [d + np.log10(m) for d in range(-6, 0) for m in range(2, 10)]
ax.set_yticks(maj); ax.set_yticks([m for m in mnr if lo < m < hi], minor=True); ax.set_yticklabels([])
ax.tick_params(axis="y", which="major", direction="in", left=True, right=True, length=47 * PT, width=3 * PT)
ax.tick_params(axis="y", which="minor", direction="in", left=True, right=True, length=22 * PT, width=3 * PT)
C = dict(SiO=(101, 194, 169), SO=(154, 211, 241), SO2=(240, 197, 102), AlO=(244, 236, 141), AlOH=(102, 171, 208), H2S=(240, 120, 111))
C = {k: tuple(c / 255 for c in v) for k, v in C.items()}
pitch = 39.2; cx = (xr - xl) / 2; x0 = cx - 3 * pitch
f = 25 / 4.0, 25 / 89.0                           # X/CO scales as 1/(12C/13C): x6.25 at 4, x0.28 at 89
vals = dict(SiO=(np.log10(2.48e-5), True), SO=(np.log10(2.29e-3), False), SO2=(np.log10(2.11e-4), False))
for i, sp in enumerate(["SiO", "SO", "SO2", "AlO", "AlOH", "H2S"]):
    xa = x0 + i * pitch; xm = xa + pitch / 2
    if sp in vals:
        v, lim = vals[sp]
        ax.bar(xa, v - lo, bottom=lo, width=pitch, align="edge", color=C[sp], edgecolor="k", linewidth=2.2, zorder=2)
        g0, g1 = max(v + np.log10(f[1]), lo), v + np.log10(f[0])
        ax.plot([xm, xm], [g0, g1], color=(0.45, 0.45, 0.45), lw=5 * PT, solid_capstyle="butt", zorder=3)
        if lim:                                   # K18-style down arrow: filled head with its tip on the bar top
            hl, hw = 34 / dec, 19.0
            ax.add_patch(Polygon([[xm - hw, v + hl], [xm + hw, v + hl], [xm, v]], closed=True, color="k", zorder=4))
    else:
        ax.text(xm, lo + 18 / dec, "n/m", rotation=90, ha="center", va="bottom", fontsize=26 * PT * 1.0, color=(0.45, 0.45, 0.45))
fig.text((xl + cx) / W, 1 - 1216 / H, "MUBLO", ha="center", va="center", fontsize=63.4 * PT)
ax.text(cx, hi - 40 / dec, "$^{12}$C/$^{13}$C = 25", ha="center", va="top", fontsize=40 * PT)
ax.text(cx, hi - 95 / dec, "gray: 4–89", ha="center", va="top", fontsize=36 * PT, color=(0.4, 0.4, 0.4))
fig.savefig(out, dpi=dpi, facecolor="white")
