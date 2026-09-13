"""Rebuild the JWST GC Treasury (GO 10678) layers for colloquium_sep2026_granada.html.

Everything is rendered onto the SAME full-bleed WCS grid the rest of the "what's
coming soon" sequence uses -- l -1.55..1.70, b +-0.688, 4"/px, 2925x1238 -- so every
PNG stacks 1:1 in the deck with no per-slide transform.

Outputs (all in assets/cmz/):
  treasury_nircam_fullbleed.png   the NIRCam mosaic itself, from the live HiPS
  treasury_miri_fullbleed.png     the MIRI F770W parallel mosaic, from the live HiPS
  fp_nircam_fullbleed.png         all 139 NIRCam pointings, coloured by visit status
  fp_miri_fullbleed.png           the same for the MIRI parallel
  fp_done_nircam_fullbleed.png    outline only, the pointings already observed
  fp_done_miri_fullbleed.png      ditto
  treasury_panel_nircam.png       the survey's starting strip, magnified, 1"/px
  treasury_panel_miri.png         ditto
  box_treasury_start.png          that strip outlined on the full-bleed grid

The panel covers a FIXED patch of sky (PANEL_L/PANEL_B below), chosen to hold both
where the survey started and the several days of pointings queued behind it, so its
placement in the deck never has to change -- re-running just fills more of it in.

Both the mosaics and the statuses are live: the HiPS grow as observations land and
the status table is STScI's own visit status, so RE-RUN THIS shortly before the talk.
It prints a provenance block; paste the counts into the slide captions if they move.

    python colloquium_assets/make_cmz_treasury_layers.py

Needs hipsrender.py (same directory) for the HiPS -> WCS resampling.
"""
import os, sys, json, warnings
import numpy as np
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hipsrender import render

import requests
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None

OUT = 'assets/cmz'
NIRCAM_HIPS = 'https://starformation.astro.ufl.edu/avm_images/jwst_gc_treasury_hips/'
MIRI_HIPS   = 'https://starformation.astro.ufl.edu/avm_images/jwst_gc_treasury_miri_hips/'
MONITOR_FP  = 'https://starformation.astro.ufl.edu/jwst-gc/monitor/footprints.json'

# HiPS order matched to the 4"/px grid: order 7 tiles are 512px => ~3.2"/px.
HIPS_ORDER = 7
# Tile orientation.  Verified by correlation, not assumed: with 'F', 99% of the
# opaque NIRCam pixels fall inside the executed NIRCam footprints; with 'A', 0%.
HIPS_MAPPING = 'F'

# full-bleed CMZ grid (inside the Ash image's coverage), recentred on l = +0.075
L0, L1 = -1.55, 1.70
B0, B1 = -0.688, 0.688
SC = 4 / 3600.

# Magnified panel.  The survey runs east to west in descending visit number: it opened
# at l = +0.70 and works down through +0.58, +0.51, +0.45 ... so this window holds the
# executed pointings and the next couple of weeks of scheduled ones.  Same projection
# and orientation as the full-bleed grid, so the panel needs no rotation in the deck.
# The MIRI parallel sits ~6' east of its NIRCam prime and reaches l = +0.852, so the
# window runs out to +0.87 to hold both instruments whole.  That is ~15 display px
# further east than the full-bleed framing reaches, so the box's left edge falls just
# off the canvas -- the same thing Sgr C's box does at the other end of the strip.
PANEL_L = (0.30, 0.87)
PANEL_B = (-0.21, 0.14)
PANEL_SC = 1 / 3600.
PANEL_ORDER = 9

# Visit status -> (display class, colour).  STScI's own vocabulary, collapsed to the
# four things an audience cares about.  Colours follow the jwst-gc monitor's sky map
# so the deck and the monitor page read the same.
STATUS_CLASS = {
    'executed': 'observed', 'archived': 'observed', 'collecting': 'observed',
    'scheduled': 'scheduled',
    'flight ready': 'planned', 'implementation': 'planned', 'planned': 'planned',
    'skipped': 'delayed', 'withdrawn': 'delayed',
}
CLASS_COLOR = {
    'observed':  (74, 222, 128),    # #4ade80
    'scheduled': (251, 191, 36),    # #fbbf24
    'delayed':   (244,  63,  94),   # #f43f5e
}
PLANNED_COLOR = {'nircam': (70, 188, 214), 'miri': (167, 139, 250)}  # #46bcd6 / #a78bfa
# Drawn last-to-first so the rare, important classes sit on top of the bulk.
DRAW_ORDER = ['planned', 'scheduled', 'observed', 'delayed']


def grid(l0=L0, l1=L1, b0=B0, b1=B1, sc=SC):
    nx = int(round((l1 - l0) / sc))
    ny = int(round((b1 - b0) / sc))
    h = fits.Header()
    for k, v in [('NAXIS', 2), ('NAXIS1', nx), ('NAXIS2', ny), ('WCSAXES', 2),
                 ('CTYPE1', 'GLON-TAN'), ('CTYPE2', 'GLAT-TAN'),
                 ('CRVAL1', (l0 + l1) / 2), ('CRVAL2', (b0 + b1) / 2),
                 ('CRPIX1', nx / 2 + 0.5), ('CRPIX2', ny / 2 + 0.5),
                 ('CDELT1', -sc), ('CDELT2', sc), ('CUNIT1', 'deg'), ('CUNIT2', 'deg')]:
        h[k] = v
    return WCS(h), nx, ny


def save(arr, name):
    p = f'{OUT}/{name}.png'
    Image.fromarray(arr).save(p, optimize=True)
    print('  %-34s %6.2f MB' % (name + '.png', os.path.getsize(p) / 1e6))


def main():
    W, NX, NY = grid()
    print('grid %dx%d  (%.1f"/px)  l %.2f..%.2f  b %+.3f..%+.3f' % (NX, NY, SC * 3600, L0, L1, B0, B1))

    def radec2px(poly):
        g = SkyCoord(ra=[p[0] for p in poly], dec=[p[1] for p in poly],
                     unit='deg', frame='icrs').galactic
        x, y = W.all_world2pix(g.l.deg, g.b.deg, 1)
        return list(zip((x - 0.5).tolist(), (NY - (y - 0.5)).tolist()))

    # ---- the mosaics themselves -------------------------------------------------
    mosaics = {}
    for key, base in [('nircam', NIRCAM_HIPS), ('miri', MIRI_HIPS)]:
        arr, nt, nu = render(base, W, (NY, NX), order=HIPS_ORDER, mapping=HIPS_MAPPING)
        mosaics[key] = arr[..., 3] > 10
        print('%-6s HiPS: %d/%d tiles present, %.4f%% of the frame carries data'
              % (key, nt, nu, 100 * mosaics[key].mean()))
        save(arr, f'treasury_{key}_fullbleed')

    # ---- footprints, coloured by visit status -----------------------------------
    fp = requests.get(MONITOR_FP, timeout=60).json()
    visits = list(fp['observed']) + list(fp['planned'])
    by_class = {c: [] for c in DRAW_ORDER}
    unknown = set()
    for v in visits:
        st = str(v.get('status', '')).strip().lower()
        cls = STATUS_CLASS.get(st)
        if cls is None:
            unknown.add(st or '(blank)')
            cls = 'planned'
        by_class[cls].append(v)
    if unknown:
        print('  !! unrecognised status(es), drawn as planned:', ', '.join(sorted(unknown)))

    print('program %s, PA_V3 %.1f deg  (%s)' % (fp['program'], fp['pa_v3'], fp['status_source']))
    print('raw STScI counts:', ', '.join('%s %d' % kv for kv in fp['status_counts'].items()))
    print('as drawn:        ', ', '.join('%s %d' % (c, len(by_class[c])) for c in DRAW_ORDER),
          ' total %d' % len(visits))

    for inst in ['nircam', 'miri']:
        im = Image.new('RGBA', (NX, NY), (0, 0, 0, 0))
        d = ImageDraw.Draw(im, 'RGBA')
        done = Image.new('RGBA', (NX, NY), (0, 0, 0, 0))
        dd = ImageDraw.Draw(done, 'RGBA')
        drawn = {}
        for cls in DRAW_ORDER:
            rgb = PLANNED_COLOR[inst] if cls == 'planned' else CLASS_COLOR[cls]
            # observed and delayed are the ones the eye must find: heavier line, more fill
            fa, lw = (95, 4) if cls in ('observed', 'delayed') else (55, 3)
            n = 0
            for v in by_class[cls]:
                for poly in (v.get(inst) or []):
                    pts = radec2px(poly)
                    if not all(np.isfinite([p[0] for p in pts])):
                        continue
                    d.polygon(pts, fill=rgb + (fa,), outline=rgb + (245,), width=lw)
                    if cls == 'observed':
                        dd.polygon(pts, outline=(74, 222, 128, 255), width=5)
                    n += 1
            drawn[cls] = n
        print('%-6s polygons drawn: %s' % (inst, drawn))
        save(np.asarray(im), f'fp_{inst}_fullbleed')
        save(np.asarray(done), f'fp_done_{inst}_fullbleed')

        # how much of what was observed has actually been reduced into the mosaic yet
        m = Image.new('L', (NX, NY), 0)
        md = ImageDraw.Draw(m)
        for v in by_class['observed']:
            for poly in (v.get(inst) or []):
                md.polygon(radec2px(poly), fill=255)
        mask = np.asarray(m) > 0
        if mask.any():
            inside = (mosaics[inst] & mask).sum() / max(mosaics[inst].sum(), 1)
            print('       mosaic vs observed footprints: %.1f%% of mosaic pixels inside, '
                  'mosaic fills %.1f%% of the observed area'
                  % (100 * inside, 100 * mosaics[inst][mask].mean()))

    panels(W, NX, NY)


def panels(W, NX, NY):
    """Magnified view of the starting strip, plus its outline on the full-bleed grid."""
    Wp, PX, PY = grid(PANEL_L[0], PANEL_L[1], PANEL_B[0], PANEL_B[1], PANEL_SC)
    print('panel %dx%d  (%.1f"/px)  l %.2f..%.2f  b %+.2f..%+.2f  aspect %.3f'
          % (PX, PY, PANEL_SC * 3600, PANEL_L[0], PANEL_L[1], PANEL_B[0], PANEL_B[1], PX / PY))
    for key, base in [('nircam', NIRCAM_HIPS), ('miri', MIRI_HIPS)]:
        arr, nt, nu = render(base, Wp, (PY, PX), order=PANEL_ORDER, mapping=HIPS_MAPPING)
        print('%-6s panel: %d/%d tiles, %.2f%% of the panel carries data'
              % (key, nt, nu, 100 * (arr[..., 3] > 10).mean()))
        save(arr, f'treasury_panel_{key}')

    # the same rectangle drawn on the wide grid, so box edges and panel edges are one edge
    im = Image.new('RGBA', (NX, NY), (0, 0, 0, 0))
    d = ImageDraw.Draw(im, 'RGBA')
    ll = [PANEL_L[1], PANEL_L[0], PANEL_L[0], PANEL_L[1]]
    bb = [PANEL_B[0], PANEL_B[0], PANEL_B[1], PANEL_B[1]]
    x, y = W.all_world2pix(np.array(ll), np.array(bb), 1)
    d.polygon(list(zip((x - 0.5).tolist(), (NY - (y - 0.5)).tolist())),
              outline=(60, 235, 60, 245), width=4)
    save(np.asarray(im), 'box_treasury_start')
    # deck placement, printed so the HTML numbers can be checked against the render
    DW, DL, DT = 2154.82, -565.41, -118.72
    s = DW / NX
    print('   box on the 1024x768 canvas: x %.1f..%.1f  y %.1f..%.1f'
          % (DL + (min(x) - 0.5) * s, DL + (max(x) - 0.5) * s,
             DT + (NY - (max(y) - 0.5)) * s, DT + (NY - (min(y) - 0.5)) * s))


if __name__ == '__main__':
    main()
