import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
import pyfits

#load utservo
raw_file = pyfits.open('/COFE/Level1/0.6/utservo_all_v0.6.fits')
ut = raw_file['GYRO_HID'].data['UT']
az = raw_file['GYRO_HID'].data['HYBRIDHEADINGANGLE']

#remove angles outside of +- pi
valid = (az < np.pi) & (az > -np.pi)
az = az[valid]
ut = ut[valid]

#just section
#sec = (ut > 32.5) & (ut < 33)
#az = az[sec]
#ut = ut[sec]
#ut_sci = ut_sci[(ut_sci > 32.5) & (ut_sci < 33)]

plt.figure()
plt.plot(ut, az, '.', label='raw')
plt.xlabel('UT')
plt.savefig('rawaz.png')

#unwrap the heading angle 
wraps = np.diff(az) < - .95 * 2 * np.pi #5% tolerance
unwrapped = az.copy()
unwrapped[1:] += np.cumsum(wraps) * np.pi * 2

typical_revlength = np.median(np.diff(ut[wraps]))

#fix single sample jumps
#second next nearer than next sample
single_sample_jumps = np.where((unwrapped[2:] - unwrapped[:-2]) < (unwrapped[1:-1] - unwrapped[:-2]))[0]+1
#create mask
continous = np.ones(len(unwrapped), dtype=np.bool)
continous[single_sample_jumps] = False
unwrapped = unwrapped[continous]
ut = ut[continous]

#fix time gaps
#all gaps longer than 1 second
h_jumps = np.diff(ut) > (5 / 3600.)
h_jumps_scaled = h_jumps.astype(np.double) 
h_jumps_scaled[h_jumps] *= np.round(np.diff(ut)[h_jumps]/typical_revlength)
unwrapped[1:] += np.cumsum(h_jumps_scaled) * np.pi * 2 

#read ut science
ut_sci_10 = pyfits.getdata('/COFE/Level1/0.6/all_10GHz_v0.6_data.fits', 'TIME')['UT']
ut_sci_15 = pyfits.getdata('/COFE/Level1/0.6/all_15GHz_v0.6_data.fits', 'TIME')['UT']

#interpolate and reset to -pi pi
fixed_az_10 = np.mod(np.interp(ut_sci_10, ut, unwrapped) + np.pi, 2*np.pi) - np.pi
fixed_az_15 = np.mod(np.interp(ut_sci_15, ut, unwrapped) + np.pi, 2*np.pi) - np.pi

plt.figure()
plt.plot(ut_sci_10, fixed_az_10, 'r.', label='fixed')
plt.xlabel('UT')
plt.savefig('fixedaz.png')

import pycfitsio as fits
fits.write('fixaz.fits', OrderedDict([
    ('10GHz', OrderedDict([('UT', ut_sci_10), ('AZ', fixed_az_10)])),
    ('15GHz', OrderedDict([('UT', ut_sci_15), ('AZ', fixed_az_15)]))
    ]))
