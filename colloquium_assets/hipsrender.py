import os, io, numpy as np, requests
from concurrent.futures import ThreadPoolExecutor
from PIL import Image
import astropy_healpix as ah
from astropy.coordinates import Galactic, SkyCoord
import astropy.units as u
Image.MAX_IMAGE_PIXELS=None
CACHE=os.path.join(os.path.dirname(os.path.abspath(__file__)),'tilecache')
os.makedirs(CACHE, exist_ok=True)
_sess=requests.Session()

def _tile(base, order, idx):
    fn=os.path.join(CACHE, f'{abs(hash(base))%10**8}_{order}_{idx}.png')
    if os.path.exists(fn):
        try: return np.asarray(Image.open(fn).convert('RGBA'))
        except Exception: pass
    d=(idx//10000)*10000
    url=f'{base.rstrip("/")}/Norder{order}/Dir{d}/Npix{idx}.png'
    try:
        r=_sess.get(url, timeout=30)
        if r.status_code!=200: return None
        open(fn,'wb').write(r.content)
        return np.asarray(Image.open(io.BytesIO(r.content)).convert('RGBA'))
    except Exception:
        return None

MAPS={
 'A': lambda dx,dy,W: (dx*W,        (1-dy)*W),
 'B': lambda dx,dy,W: (dy*W,        (1-dx)*W),
 'C': lambda dx,dy,W: ((1-dx)*W,    dy*W),
 'D': lambda dx,dy,W: ((1-dy)*W,    dx*W),
 'E': lambda dx,dy,W: (dx*W,        dy*W),
 'F': lambda dx,dy,W: (dy*W,        dx*W),
 'G': lambda dx,dy,W: ((1-dx)*W,    (1-dy)*W),
 'H': lambda dx,dy,W: ((1-dy)*W,    (1-dx)*W),
}

def render(base, wcs, shape, order, tile_w=512, mapping='A', frame='galactic', workers=12):
    NY,NX=shape
    jj,ii=np.mgrid[0:NY,0:NX]
    lon,lat=wcs.all_pix2world(ii+1.0, NY-jj, 1)      # image row -> FITS y
    if frame=='galactic':
        L,B=lon,lat
    else:
        c=SkyCoord(ra=lon,dec=lat,unit='deg',frame='icrs').galactic; L,B=c.l.deg,c.b.deg
    hp=ah.HEALPix(nside=2**order, order='nested', frame=Galactic())
    idx,dx,dy=hp.lonlat_to_healpix(L*u.deg, B*u.deg, return_offsets=True)
    col,row=MAPS[mapping](dx,dy,tile_w)
    col=np.clip(col.astype(int),0,tile_w-1); row=np.clip(row.astype(int),0,tile_w-1)
    out=np.zeros((NY,NX,4),np.uint8)
    uniq=np.unique(idx)
    tiles={}
    with ThreadPoolExecutor(workers) as ex:
        for t,arr in zip(uniq, ex.map(lambda t:_tile(base,order,int(t)), uniq)):
            if arr is not None: tiles[int(t)]=arr
    for t,arr in tiles.items():
        m=idx==t
        out[m]=arr[row[m],col[m]]
    return out, len(tiles), len(uniq)
