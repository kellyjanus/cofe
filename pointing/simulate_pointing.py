import datetime
import healpy as hp
from collections import OrderedDict
import numpy as np
from scipy import constants
from pointingtools import compute_parallactic_angle, altaz2ha 

import pycfitsio as fits
import ephem

NSIDE = 128
NPIX = hp.nside2npix(NSIDE)

lat = np.radians([-23])
lon = np.radians([0])
alt = [30000]

def conv(i, azimuth, elevation, utc):
    observer = ephem.Observer()
    observer.lon = lon[i]
    observer.lat = lat[i]
    observer.elevation = alt[i]  
    observer.date = utc
    return observer.radec_of(azimuth, elevation)

START_DATE = datetime.datetime(2012, 8, 1)

rotation_speed = np.radians(-1 * 360/60)
ut = np.arange(0., 48., 1/(30.*3600.))
az = rotation_speed * (ut * 3600.) % (2*np.pi)

ra = []
dec = []
el = np.radians(45)
for i in range(0, len(ut)):
    ra_i, dec_i = conv(0, az[i], el, START_DATE+datetime.timedelta(ut[i]/24.))
    ra.append(ra_i)
    dec.append(dec_i)

ra = np.array(ra)
dec = np.array(dec)
ha = altaz2ha(el, az, lat)
psi = compute_parallactic_angle(ha, lat, dec)

pix = hp.ang2pix(NSIDE, np.pi/2-dec, ra, nest=False)

from mapmaking import pix2map
hit = hp.ma(pix2map(pix, NSIDE))
hit.mask = hit == 0
hp.mollview(np.log10(hit).filled(),min=0,max=3,unit='log10 hits',title='Hitmap 48hrs Alice Springs')
