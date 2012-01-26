import datetime
import numpy as np
from scipy import constants

import pycfitsio as fits
import ephem

def freq2wavelength(freq):
        """Freq [GHz] to wavelength [microns]"""
        return constants.c / freq / 1e3

freq = 15
pointing_file = fits.open('data/v5_all_%dghz_pointing.fits' % freq)
data_file = fits.open('data/all_%dGHz_v0.5_data.fits' % freq)
servo_file = fits.open('data/utservo_all_v0.5.fits')

# pointing channel is the column in the pointing file
pnt_ch = 1

MISSIONSTART = 16.8 #from altitude
MISSIONEND = 36.84 #from issue with latitude
#MISSIONEND = 16.8 + 1./60
#first sun xsing
#MISSIONSTART = 17 + 9/60.
#MISSIONEND = 17 + 13/60.
START_DATE = datetime.date(2011, 9, 17)
START_JULIAN = ephem.julian_date(START_DATE) + MISSIONSTART/24.


#azimuth/elevation
ut = data_file['TIME'].read_column('UT') 
good = (ut > MISSIONSTART) & (ut < MISSIONEND)
ut = ut[good]
az = np.radians(pointing_file[0].read_column('az%d' % pnt_ch)[good])
el = np.radians(pointing_file[0].read_column('el%d' % pnt_ch)[good])

#sanitize altitude, latitude and longitude
servo_ut = servo_file['GYRO_HID'].read_column('UT')
TOL = 1./60.
servo_good = (servo_ut > MISSIONSTART - TOL) & (servo_ut < MISSIONEND + TOL)
servo_alt = servo_file['GYRO_HID'].read_column('HYBRIDALTITUDE')
good_alt = servo_good & (servo_alt > 4000) & (servo_alt < 1e5)
alt = np.interp(ut, servo_ut[good_alt], servo_alt[good_alt])
servo_lat = np.degrees(servo_file['GYRO_HID'].read_column('HYBRIDLATITUDE'))
good_lat =  servo_good &(servo_lat < 37) & (servo_lat > 32)
lat = np.radians(np.interp(ut, servo_ut[good_lat], servo_lat[good_lat]))
servo_lon = np.degrees(servo_file['GYRO_HID'].read_column('HYBRIDLONGITUDE'))
good_lon = servo_good & (servo_lon > -114) & (servo_lon < -86)
lon = np.radians(np.interp(ut, servo_ut[good_lon], servo_lon[good_lon]))

J0 = ephem.julian_date(0)
ra = []
dec = []
ut_out = []
import time
#def conv(i):
#    ra, dec = tpmu.convert(x=az[i], y=el[i], s1=19, s2=7, epoch=2451545.0,
#                  equinox=2451545.0, timetag=START_JULIAN + (ut[i]-MISSIONSTART)/24., 
#                  lon=lon[i], lat=lat[i], alt=alt[i],
#                  T=273.15, P=1023.25,
#                  H=0.0, W=freq2wavelength(freq))
#    return ra,dec

def conv(i):

    observer = ephem.Observer()
    observer.lon = lon[i]
    observer.lat = lat[i]
    observer.elevation = alt[i]  
    observer.date = START_JULIAN + (ut[i]-MISSIONSTART)/24. - J0

    return observer.radec_of(az[i], el[i])

#def conv(i):
#    v6 = convert.cat2v6(alpha = az[i], delta = el[i])
#    start_clock = time.clock()
#    v6c = convert.convertv6(v6=v6,
#        utc=START_JULIAN + (ut[i]-MISSIONSTART)/24.,# delta_at=-999, delta_ut=-999,
#        s1=tpm.TPM_S19, s2=tpm.TPM_S07,
#        epoch=tpm.J2000, equinox=tpm.J2000,
#        lon=lon[i], lat=lat[i], alt=alt[i],
#        xpole=0.0, ypole=0.0,
#        T=273.15, P=1013.25, H=0.0, wavelength=freq2wavelength(freq))
#    print("TIME[%d]:%.2g s" % (i, time.clock() - start_clock))
#    cat = convert.v62cat(v6c)
#    return cat['alpha'], cat['delta']
#for i in range(0, 10):

for i in range(0, len(ut)):
    print(i * 100./len(ut))
    ut_out.append(ut[i])
    ra_i, dec_i = conv(i)
    ra.append(ra_i)
    dec.append(dec_i)

np.save('ut_out_%d_ch%d' % (freq, pnt_ch), np.array(ut_out))
np.save('ra_%d_ch%d' % (freq, pnt_ch)  ,np.array(ra) )
np.save( 'dec_%d_ch%d' % (freq, pnt_ch),np.array(dec))
