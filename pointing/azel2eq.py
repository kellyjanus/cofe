import datetime
from collections import OrderedDict
import numpy as np
from scipy import constants
from pointingtools import compute_parallactic_angle, altaz2ha 

import pycfitsio as fits
import ephem

def freq2wavelength(freq):
        """Freq [GHz] to wavelength [microns]"""
        return constants.c / freq / 1e3

freq = 10
mag = True

pointing_filename = 'data/all_%dghz_pointing.fits' % freq
if mag:
    pointing_filename = pointing_filename.replace('.fits','_mag.fits')
pointing_file = fits.open(pointing_filename)

channels = list(set([int(c.translate(None, 'AZEL')) for c in pointing_file[0].dtype.names if c != 'UT']))
data_file = fits.open('data/all_%dGHz_data_cal.fits' % freq)
servo_file = fits.open('data/utservo.fits')

# pointing channel is the column in the pointing file

MISSIONSTART = 16.8 #from altitude
MISSIONEND = 36.76 #from issue with latitude
#MISSIONEND = 16.8 + 1./60
#first sun xsing
#MISSIONSTART = 17 + 9/60.
#MISSIONEND = 17 + 13/60.
START_DATE = datetime.datetime(2011, 9, 17)

#azimuth/elevation
ut = data_file['TIME'].read_column('UT') 
good = (ut > MISSIONSTART) & (ut < MISSIONEND)
ut = ut[good]

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

def conv(i, azimuth, elevation, utc):
    observer = ephem.Observer()
    observer.lon = lon[i]
    observer.lat = lat[i]
    observer.elevation = alt[i]  
    observer.date = utc
    return observer.radec_of(azimuth, elevation)

out_filename = 'data/eq_pointing_%d.fits' % (freq)
if mag:
    out_filename = out_filename.replace('.fits','_mag.fits')
with fits.create(out_filename) as f:
    f.write_HDU("TIME", OrderedDict({'UT': ut}))
    for pnt_ch in channels:
        print("Channel %d" % pnt_ch)
        az = np.radians(pointing_file[0].read_column('az%d' % pnt_ch)[good])
        el = np.radians(pointing_file[0].read_column('el%d' % pnt_ch)[good])

        ra = []
        dec = []
        for i in range(0, len(ut)):
            if i % 100000 == 0:
                print("%.1d perc" % (i * 100./len(ut)))
            ra_i, dec_i = conv(i, az[i], el[i], START_DATE+datetime.timedelta(ut[i]/24.))
            ra.append(ra_i)
            dec.append(dec_i)

        ra = np.array(ra)
        dec = np.array(dec)
        ha = altaz2ha(el, az, lat)
        psi = compute_parallactic_angle(ha, lat, dec)

        pnt = OrderedDict()
        pnt['THETA'] = np.pi/2 - dec
        pnt['PHI'] = ra
        pnt['PSI'] = psi
        f.write_HDU("CHANNEL_%d" % pnt_ch, pnt)
