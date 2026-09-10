#!/usr/bin/env python3
"""
Sigma_SFR vs Sigma_gas: the local single-cloud star formation relation, with CMZ clouds on top.

Reproduces Pokhrel et al. (2021, ApJL 912, L19) Fig. 2a and overlays the Galactic Centre
clouds we have numbers for.  Rerun this to rebuild the figure as new data arrive -- edit
CMZ_CLOUDS below and nothing else.

    python3 colloquium_assets/make_sf_relation_figure.py

Outputs colloquium_assets/sf_relation_cmz.png (and .pdf).

The local cloud tracks come from colloquium_assets/pokhrel2021_tracks.csv, digitised from
Fig. 2a of the paper (~/Dropbox (Personal)/Papers/Pokhrel2021_UniformEpsFF.pdf) by colour-
keying each curve against its legend swatch.  Ten of the twelve clouds are recoverable that
way; Orion-A (black) and Mon R2 (grey) are drawn in colours indistinguishable from the axes
and the fit line, so they are not included.  Rerun with --digitise to rebuild the CSV.

UNITS THROUGHOUT: Sigma_gas [M_sun / pc^2],  Sigma_SFR [M_sun / pc^2 / Myr].
"""
import argparse, csv, os, sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
TRACKS = os.path.join(HERE, "pokhrel2021_tracks.csv")
SGRB2 = os.path.join(HERE, "ginsburg2018_fig17a.csv")
PAPER = os.path.expanduser("~/Dropbox (Personal)/Papers/Pokhrel2021_UniformEpsFF.pdf")
OUT = os.path.join(HERE, "sf_relation_cmz")

# ----------------------------------------------------------------------------------------
# Pokhrel et al. 2021.  Table 1 and Section 4.1.
#   log Sigma_SFR = a log Sigma_gas + b ; sample-average a = 2.00 +/- 0.27, b = -4.11 +/- 0.80.
#   Cloud-to-cloud scatter sigma = 0.30 dex, evaluated at Sigma_gas = 10^2.5 (Table 1 note).
#   Fig. 2 caption: dark band = +/-1 sigma, light band = +/-2 sigma.
#   Sigma_SFR = eps_ff Sigma_gas / t_ff with median eps_ff = 0.026 (Section 4.2).
#   Their YSO -> SFR recipe (Section 3): SFR = N_PS * 0.5 Msun / 0.5 Myr, i.e. protostellar
#   (Class 0/I) counts only.
# THE FIT RANGE MATTERS: Section 4.1 fits only below N(H2) = 3e22 cm^-2, i.e.
#   Sigma_gas <~ 676 Msun/pc^2.  Everything from the CMZ is beyond that.
# ----------------------------------------------------------------------------------------
POKHREL_SLOPE, POKHREL_INTERCEPT = 2.00, -4.11
POKHREL_SCATTER_DEX = 0.30
POKHREL_FIT_MAX_SIGMA = 676.0
POKHREL_EPS_FF = 0.026

# Conversions used for the CMZ points, kept here so they are easy to see and change.
PROTOSTELLAR_LIFETIME_MYR = 0.5     # Pokhrel's t_PS, reused for the CMZ so the axes mean the same thing

# Sgr B2 is drawn as the whole Fig. 17a point cloud, not one point: Sigma_SFR = Sigma_* / age.
# Ginsburg+2018 use 0.74 Myr (time since pericentre).  0.5 Myr is used here because it matches
# Pokhrel's t_PS, so both axes mean the same thing, AND it is the more conservative choice --
# a shorter age raises Sigma_SFR and shrinks the gap.  0.1 Myr is the extreme: even then the
# cloud stays below the extrapolated relation (see the printout).
SGRB2_AGE_MYR = 0.5
CORE_TO_STAR_EFFICIENCY = 0.3       # Ginsburg+2018 Eq. 3 "star formation efficiency of a core"

# N(H2) [cm^-2] -> Sigma_gas [Msun/pc^2] with mu = 2.8 amu (the convention Ginsburg, Barnes and
# Walker all use).  Pokhrel's 2 m_H N(H2)/X with X = 0.71 agrees with this to 1%.
N_H2_TO_SIGMA = 2.2419e-20

# ----------------------------------------------------------------------------------------
# CMZ clouds.  ADD ROWS HERE AS DATA ARRIVE.
#   sigma_gas / sigma_sfr are linear (not log).  Set limit=True for an upper limit.
#   sigma_sfr_lo/hi give a range to draw as an error bar (None to skip).
#   Every row must carry its own derivation in `note` -- if you cannot write the derivation,
#   the point does not belong on the plot.
# ----------------------------------------------------------------------------------------
CMZ_CLOUDS = [
    dict(
        name="The Brick",
        sigma_gas=5.1e3, sigma_gas_lo=1.0e23 * 2.2419e-20, sigma_gas_hi=6.4e23 * 2.2419e-20,
        sigma_sfr=10.6, sigma_sfr_lo=None, sigma_sfr_hi=None,
        limit=True, colour="#1f77b4", marker="v", size=230, label_offset=(-16, -34), label_ha="right",
        ref="Walker+ 2021",
        note=("CLOUD-SCALE aperture, to match how Pokhrel measure (all the protostars in a cloud "
              "against all the gas above a contour) and so that the Brick and Cloud E/F are "
              "comparable to each other.  "
              "Sigma_gas: >1e5 Msun within R = 2-3 pc (Sec 1, quoting Immer+2012, Longmore+2012, "
              "Walker+2015); R = 2.5 pc -> A = 19.6 pc^2 -> 5.1e3 Msun/pc^2.  "
              "Sigma_SFR: DERIVED -- Walker+ quote no SFR and no limit.  Generous ceiling on the "
              "STELLAR content: every Table 3 core at its 22 K mass, source 1 counted whole "
              "(64.2 Msun) rather than split, = 104 Msun of core GAS, times a core-to-star "
              "efficiency of 1.0 -- i.e. every gram already a star.  / 0.5 Myr / 19.6 pc^2 "
              "-> 10.6 Msun/pc^2/Myr.  "
              "Horizontal bar: the range of column the cloud actually covers, from the CMZoom "
              "threshold N(H2) = 1e23 cm^-2 (Walker+ 2018 Sec 2.1: 'designed to target all regions "
              "within the CMZ that lie above a column density threshold of ~1e23 cm^-2') up to the "
              "densest structure in the cloud, the maser core -- 72 Msun in R = 0.04 pc "
              "(Rathborne+ 2015 via Walker+ 2021 Sec 4.2) = 1.4e4 Msun/pc^2 = 6.4e23 cm^-2.  "
              "These cores are genuinely protostellar: 9 of 18 drive SiO outflows.  "
              "If instead you use the 19-arcsec ALMA field alone (7e3 Msun, 1.75 pc^2) the point "
              "moves to (4.0e3, 119) -- same cloud, aperture ten times smaller.  Fiducial "
              "efficiency 0.3 on the 44 Msun low-temperature core sum gives 1.3 here.  "
              "CROSS-CHECK (Tang+ 2021, MNRAS 505, 2392, Fig. 10 -- one uniform AzTEC+Herschel "
              "dust-SED N(H2) map of the whole CMZ, so the Brick and Cloud E/F are measured the "
              "same way): reading the map in matched R = 2.4 pc apertures gives mean lg N(H2) = "
              "22.90 for the Brick against 23.03/23.01 for clouds e/f.  Cloud E/F is the denser "
              "of the two by ~0.13 dex, which is the same ordering the Walker+ 2018 masses give "
              "(5.1e3 vs 6.0e3), and by a slightly larger margin -- so the separation plotted "
              "here is the conservative one."),
    ),
    dict(
        name="Cloud E/F",
        sigma_gas=6.03e3, sigma_gas_lo=1.0e23 * 2.2419e-20, sigma_gas_hi=3.7e24 * 2.2419e-20,
        sigma_sfr=1.96, sigma_sfr_lo=None, sigma_sfr_hi=None,
        limit=True, colour="#2ca02c", marker="v", size=230, label_offset=(18, 22),
        ref="Barnes+ 2019, Walker+ 2018",
        note=("Sigma_gas: from Walker+ 2018 Table 1, which tabulates the dust ridge clouds on a "
              "consistent footing (values from Walker+ 2015).  Cloud e: 11.2e4 Msun, R = 2.4 pc; "
              "cloud f: 7.3e4 Msun, R = 2.0 pc.  Barnes+ treat e and f as one cloud, so combined: "
              "1.85e5 Msun over 30.7 pc^2 -> 6.0e3 Msun/pc^2.  "
              "The Immer+ 2012 column of the same table (15.3e4/4.5 pc, 7.2e4/2.7 pc) gives "
              "2.6e3 instead, so treat the horizontal position as good to a factor of ~2.  "
              "NOT used: Barnes' intro prose '~1e5 Msun, radii of ~1 pc'.  Multiplied together "
              "those give 3.2e4 Msun/pc^2, five to thirteen times too high -- they are "
              "order-of-magnitude statements, and Barnes never tabulate a cloud-scale mass or "
              "radius.  That error put Cloud E/F to the right of most of Sgr B2, which is wrong: "
              "it should sit at lower column density than most of the Sgr B2 cloud.  "
              "Sigma_SFR: ceiling set by what is OBSERVED.  No HII regions, no protostellar "
              "source catalogue; the only star formation signature in either cloud is an H2O plus "
              "Class II CH3OH maser pair at one spot in the south.  A Class II methanol maser marks "
              "a high-mass protostar, so allow one such star plus companions, <= 30 Msun -- "
              "generous, and bounded above by the Brick's ENTIRE counted protostellar content "
              "(104 Msun in 18 cores), since E/F shows strictly less.  / 0.5 Myr / 30.7 pc^2 "
              "-> 2.0 Msun/pc^2/Myr.  "
              "Horizontal bar: from the same 1e23 cm^-2 threshold up to the peak column Barnes+ "
              "measure, 3.7e24 cm^-2 towards the south of Cloud E/F (Sec 3.1) = 8.3e4 Msun/pc^2.  "
              "The bar is drawn flat, at the cloud-averaged limit.  A contour-by-contour treatment "
              "would tilt it upward to the right (smaller area at higher contour), the way Pokhrel's "
              "clouds are curves rather than points -- so the flat bar is the conservative rendering "
              "at the high-column end.  "
              "NOT used: Barnes' 'core region could be up to ~600 Msun ... capable of forming one "
              "or several high-mass stars (~25 per cent)'.  That 150 Msun is star formation "
              "potential, not stars that exist -- these cores are starless.  "
              "See the Brick entry for the Tang+ 2021 cross-check that puts Cloud E/F above the "
              "Brick in column density on a single uniform map."),
    ),
    # dict(name="Clouds C & D", ... ),   # Gramze+ in prep
    # dict(name="Sgr C", ... ),          # Crowe+ 2023, Lu+ 2021
]


def digitise(csv_path=TRACKS, pdf=PAPER):
    """Rebuild pokhrel2021_tracks.csv from the paper by colour-keying Fig. 2a."""
    import fitz
    from PIL import Image
    doc = fitz.open(pdf)
    xref = None
    for page in doc:
        for img in page.get_images(full=True):
            pix = fitz.Pixmap(doc, img[0])
            if pix.width > 1800 and 600 < pix.height < 800:
                xref = img[0]
                break
        if xref:
            break
    if xref is None:
        sys.exit("could not find the two-panel figure in %s" % pdf)
    pix = fitz.Pixmap(doc, xref)
    if pix.n > 4:
        pix = fitz.Pixmap(fitz.csRGB, pix)
    tmp = os.path.join(HERE, "_pokhrel_fig.png")
    pix.save(tmp)
    a = np.asarray(Image.open(tmp).convert("RGB")).astype(int)
    os.remove(tmp)
    # panel (a) axes box, measured from the long ink columns/rows of the extracted image
    X0, X1, Y0, Y1 = 119, 970, 3, 597
    XL, XR, YT, YB = 1.8, 3.4, 2.5, -0.5
    clouds = [("Ophiuchus", (57, 236, 209)), ("Perseus", (47, 169, 244)),
              ("Orion-B", (233, 47, 47)), ("Aquila-North", (202, 248, 50)),
              ("Aquila-South", (185, 118, 72)), ("NGC 2264", (252, 114, 106)),
              ("S140", (3, 68, 3)), ("AFGL 490", (236, 12, 217)),
              ("Cep OB3", (251, 205, 53)), ("Cygnus-X", (203, 167, 240))]
    sub = a[Y0:Y1 + 1, X0:X1 + 1]
    satmap = sub.max(2) - sub.min(2)
    rows = []
    for name, c in clouds:
        d = np.sqrt(((sub - np.array(c)) ** 2).sum(2))
        m = (d < 45) & (satmap > 35)          # colour match AND saturated: keeps out axes/bands
        for j in range(sub.shape[1]):
            ys = np.where(m[:, j])[0]
            if len(ys) == 0:
                continue
            y = np.median(ys)
            rows.append((name,
                         round(XL + j / (X1 - X0) * (XR - XL), 4),
                         round(YT + y / (Y1 - Y0) * (YB - YT), 4)))
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cloud", "log_sigma_gas", "log_sigma_sfr"])
        w.writerows(rows)
    print("wrote %s (%d points, %d clouds)" % (csv_path, len(rows), len(clouds)))


def digitise_sgrb2(csv_path=SGRB2,
                   pdf=os.path.expanduser("~/Dropbox (Personal)/Papers/Ginsburg2018_Distributed.pdf")):
    """Rebuild ginsburg2018_fig17a.csv: the Fig. 17a point cloud, Sigma_* vs Sigma_gas.

    The embedded raster spans panel (a)'s left edge to panel (b)'s right edge and the full axes
    height.  Panel (a) is its first 1408 columns.  Axis calibration (log Sigma_gas 3..5,
    log Sigma_* 0..5) was verified against the figure's own Gutermuth alpha=1 model lines: at
    Sigma_gas = 1e4 the t = 0.74 and 0.01 Myr lines are predicted at page rows 678 and 989 and
    found at 677 and 993.
    """
    import fitz
    from PIL import Image
    from scipy import ndimage
    doc = fitz.open(pdf)
    page = doc[18]
    pix = fitz.Pixmap(doc, page.get_images(full=True)[0][0])
    if pix.n > 4:
        pix = fitz.Pixmap(fitz.csRGB, pix)
    tmp = os.path.join(HERE, "_g18_fig17.png")
    pix.save(tmp)
    a = np.asarray(Image.open(tmp).convert("RGB")).astype(int)
    os.remove(tmp)
    H = a.shape[0]
    PW = 1408
    sub = a[:, :PW]
    dark = (sub.max(2) < 110) & ((sub.max(2) - sub.min(2)) < 45)   # black markers only
    lab, n = ndimage.label(dark)
    sizes = ndimage.sum(dark, lab, range(1, n + 1))
    keep = [i + 1 for i, sz in enumerate(sizes) if 6 <= sz <= 400]
    cen = np.array(ndimage.center_of_mass(dark, lab, keep))
    y, x = cen[:, 0], cen[:, 1]
    m = x < PW * 0.965        # drop the saturated-column lower-limit triangles at the right edge
    lg = 3.0 + x[m] / PW * 2.0
    lst = 5.0 - y[m] / H * 5.0
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["log_sigma_gas", "log_sigma_star"])
        for i in range(len(lg)):
            w.writerow([round(lg[i], 4), round(lst[i], 4)])
    print("wrote %s (%d points)" % (csv_path, len(lg)))


def load_sgrb2(csv_path=SGRB2):
    d = np.genfromtxt(csv_path, delimiter=",", names=True)
    return d["log_sigma_gas"], d["log_sigma_star"]


def load_tracks(csv_path=TRACKS):
    d = defaultdict(list)
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            d[r["cloud"]].append((float(r["log_sigma_gas"]), float(r["log_sigma_sfr"])))
    return {k: np.array(sorted(v)) for k, v in d.items()}


def make_figure(out=OUT, xlim=(1.8, 5.0), ylim=(-0.6, 4.0)):
    tracks = load_tracks()
    fig, ax = plt.subplots(figsize=(9.2, 6.9), dpi=150)

    x = np.linspace(xlim[0], xlim[1], 400)
    y = POKHREL_SLOPE * x + POKHREL_INTERCEPT
    ax.fill_between(x, y - 2 * POKHREL_SCATTER_DEX, y + 2 * POKHREL_SCATTER_DEX,
                    color="0.85", zorder=0)
    ax.fill_between(x, y - POKHREL_SCATTER_DEX, y + POKHREL_SCATTER_DEX,
                    color="0.70", zorder=0)
    ax.plot(x, y, "k--", lw=2.0, zorder=4)

    # where the relation was actually fitted, versus where we are extrapolating it
    xfit = np.log10(POKHREL_FIT_MAX_SIGMA)
    ax.axvspan(xfit, xlim[1], color="w", alpha=0.45, zorder=1)
    ax.axvline(xfit, color="0.45", lw=1.1, ls=":", zorder=5)
    ax.text(xfit - 0.05, ylim[0] + 0.10, "Pokhrel+ fit range", ha="right", va="bottom",
            fontsize=11, color="0.35")
    ax.text(xfit + 0.05, ylim[0] + 0.10, "extrapolated", ha="left", va="bottom",
            fontsize=11, color="0.35")

    for name, v in sorted(tracks.items()):
        ax.plot(v[:, 0], v[:, 1], lw=1.5, alpha=0.9, zorder=3)

    # Sgr B2: the whole Fig. 17a cloud, shifted down by the assumed age
    sg, sstar = load_sgrb2()
    ssfr = sstar - np.log10(SGRB2_AGE_MYR)
    from scipy.stats import gaussian_kde
    k = gaussian_kde(np.vstack([sg, ssfr]))
    gx_, gy_ = np.mgrid[sg.min() - .3:sg.max() + .3:120j, ssfr.min() - .3:ssfr.max() + .3:120j]
    z = k(np.vstack([gx_.ravel(), gy_.ravel()])).reshape(gx_.shape)
    zs = np.sort(z.ravel())[::-1]
    cum = np.cumsum(zs) / zs.sum()
    lv = sorted(zs[np.searchsorted(cum, f)] for f in (0.90, 0.50))
    ax.contourf(gx_, gy_, z, levels=lv + [z.max()], colors=["#d62728", "#b01d1e"], alpha=0.20, zorder=6)
    ax.contour(gx_, gy_, z, levels=lv, colors="#8c1416", linewidths=1.6, zorder=6)
    ax.scatter(sg, ssfr, s=3, c="#8c1416", alpha=0.35, lw=0, zorder=6)
    # cloud-averaged Sgr B2: the SAME kind of measurement as the two limits, so the comparison
    # at fixed Sigma_gas is like for like (the blob above is resolved on 0.25 pc cells).
    ax.scatter([np.log10(6667.0)], [np.log10(216.0)], s=300, marker="*",
               facecolors="none", edgecolors="#8c1416", linewidths=2.2, zorder=9)
    ax.plot([np.log10(6667.0)] * 2, [np.log10(156.0), np.log10(276.0)],
            color="#8c1416", lw=2.4, zorder=9)
    ax.annotate("Sgr B2, cloud-averaged", (np.log10(6667.0) + 0.06, np.log10(216.0)),
                ha="left", va="center", fontsize=11, color="#8c1416", zorder=9,
                bbox=dict(facecolor="white", alpha=0.80, edgecolor="none", pad=1.6))
    ax.annotate("Sgr B2\nGinsburg+ 2018", #(Fig. 17a / %.1f Myr)" % SGRB2_AGE_MYR,
                (4.30, 3.30), ha="center", fontsize=13, color="#a11d1f", fontweight="bold",
                zorder=9, bbox=dict(facecolor="white", alpha=0.80, edgecolor="none", pad=2.0))

    for c in CMZ_CLOUDS:
        gx, gy = np.log10(c["sigma_gas"]), np.log10(c["sigma_sfr"])
        ax.scatter([gx], [gy], s=c["size"], marker=c["marker"], c=c["colour"],
                   edgecolors="k", linewidths=1.2, zorder=8)
        if c.get("sigma_gas_lo"):
            ax.plot([np.log10(c["sigma_gas_lo"]), np.log10(c["sigma_gas_hi"])], [gy, gy],
                    color=c["colour"], lw=3.4, solid_capstyle="butt", zorder=7)
            for e in ("sigma_gas_lo", "sigma_gas_hi"):
                ax.plot([np.log10(c[e])] * 2, [gy - 0.07, gy + 0.07],
                        color=c["colour"], lw=3.4, zorder=7)
        if c["limit"]:
            ax.annotate("", xy=(gx, gy - 0.42), xytext=(gx, gy - 0.06),
                        arrowprops=dict(arrowstyle="-|>", color=c["colour"], lw=2.4), zorder=8)
        elif c.get("sigma_sfr_lo"):
            ax.plot([gx, gx], [np.log10(c["sigma_sfr_lo"]), np.log10(c["sigma_sfr_hi"])],
                    color=c["colour"], lw=2.6, solid_capstyle="butt", zorder=7)
        ax.annotate("%s\n%s" % (c["name"], c["ref"]), (gx, gy),
                    textcoords="offset points", xytext=c.get("label_offset", (16, 10)),
                    ha=c.get("label_ha", "left"),
                    fontsize=13, color=c["colour"], fontweight="bold", zorder=9,
                    bbox=dict(facecolor="white", alpha=0.80, edgecolor="none", pad=2.0))

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$\log\ \Sigma_{\rm gas}\ (M_\odot\ {\rm pc}^{-2})$", fontsize=15)
    ax.set_ylabel(r"$\log\ \Sigma_{\rm SFR}\ (M_\odot\ {\rm pc}^{-2}\ {\rm Myr}^{-1})$", fontsize=15)
    ax.tick_params(labelsize=13)

    handles = [Line2D([], [], color="0.3", ls="--", lw=2,
                      label=r"Pokhrel+ 2021: $\Sigma_{\rm SFR}\propto\Sigma_{\rm gas}^{2.0}$"),
               Line2D([], [], color="0.70", lw=9, label=r"$\pm1\sigma$, $\pm2\sigma$ cloud-to-cloud"),
               Line2D([], [], color="C0", lw=1.6, label="local clouds (10 of 12)")]
    leg = ax.legend(handles=handles, loc="upper left", fontsize=12, framealpha=0.92)
    leg.set_zorder(10)

    fig.tight_layout()
    fig.savefig(out + ".png")
    fig.savefig(out + ".pdf")
    print("wrote %s.png / .pdf" % out)

    print("\nprovenance")
    print("-" * 78)
    med_g, med_s = np.median(sg), np.median(ssfr)
    print("Sgr B2       Fig. 17a cloud, %d points, Sigma_* / %.2f Myr" % (len(sg), SGRB2_AGE_MYR))
    print("    median (log Sgas, log Ssfr) = (%.2f, %.2f); relation predicts %.2f -> %.2f dex below"
          % (med_g, med_s, POKHREL_SLOPE * med_g + POKHREL_INTERCEPT,
             POKHREL_SLOPE * med_g + POKHREL_INTERCEPT - med_s))
    print("    cloud-averaged, for comparison with the two limits: 1.5e6 Msun and "
          "SFR 0.035-0.062 Msun/yr over the 15x15 pc field -> (3.82, 2.19-2.44)")
    for age in (1.0, 0.5, 0.1):
        m = np.median(sstar - np.log10(age))
        print("      age %.2f Myr -> median log Ssfr %.2f, %.2f dex below the relation"
              % (age, m, POKHREL_SLOPE * med_g + POKHREL_INTERCEPT - m))
    print()
    for c in CMZ_CLOUDS:
        lg, ls = np.log10(c["sigma_gas"]), np.log10(c["sigma_sfr"])
        pred = POKHREL_SLOPE * lg + POKHREL_INTERCEPT
        print("%-12s log Sgas=%.2f  log Ssfr=%.2f%s  [%s]"
              % (c["name"], lg, ls, "  (upper limit)" if c["limit"] else "", c["ref"]))
        print("    relation predicts %.2f -> %.2f dex below (2 sigma = %.2f dex)"
              % (pred, pred - ls, 2 * POKHREL_SCATTER_DEX))
        for line in c["note"].split("  "):
            if line.strip():
                print("    " + line.strip())
        print()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--digitise", action="store_true",
                   help="rebuild pokhrel2021_tracks.csv from the paper PDF")
    args = p.parse_args()
    if args.digitise or not os.path.exists(TRACKS):
        digitise()
    if args.digitise or not os.path.exists(SGRB2):
        digitise_sgrb2()
    make_figure()
