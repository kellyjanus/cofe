import datetime
import time
import numpy as np

from pytpm import convert, tpm 

az = 3.30084818 #rad
el = 0.94610742 #rad
lat = 34.64 #deg
lon = -103.7 #deg
alt = 35800.26 #m
ut = 2455822.20000367 #julian date

for i in range(0, 1):
    v6 = convert.cat2v6(alpha = az, delta = el, pma=0.0, pmd=0.0, px=0.0, rv=0, C=tpm.CJ)
    start_clock = time.clock()
    v6c = convert.convertv6(v6=v6,
        utc=ut,
        s1=tpm.TPM_S19, s2=tpm.TPM_S07,
        epoch=tpm.J2000, equinox=tpm.J2000,
        lon=lon, lat=lat, alt=alt,
        xpole=0.0, ypole=0.0,
        T=273.15, P=1013.25, H=0.0, wavelength=19986.16386)
    print("TIME[%d]:%.2g s" % (i, time.clock() - start_clock))
    cat = convert.v62cat(v6c)
    print(np.degrees([cat['alpha'], cat['delta']]))
