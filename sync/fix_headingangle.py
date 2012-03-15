import numpy as np
from smooth import smooth
from collections import OrderedDict
#import matplotlib.pyplot as plt
import pyfits
import os

#load utservo
folder = os.getcwd()
utservo_filename = os.path.join(folder,'utservo.fits')
print utservo_filename
raw_file = pyfits.open(utservo_filename)
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

#plt.figure()
#plt.plot(ut, az, '.', label='raw')
#plt.xlabel('UT')
#plt.savefig('rawaz.png')

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

#smooth
#unwrapped = smooth(unwrapped, 30)

#read ut science
ut_sci_10 = pyfits.getdata(os.path.join(folder, 'all_10GHz_data.fits'), 'TIME')['UT']
ut_sci_15 = pyfits.getdata(os.path.join(folder, 'all_15GHz_data.fits'), 'TIME')['UT']

#interpolate and reset to -pi pi
fixed_az_10 = np.mod(smooth(np.interp(ut_sci_10, ut, unwrapped),30) + np.pi, 2*np.pi) - np.pi
fixed_az_15 = np.mod(smooth(np.interp(ut_sci_15, ut, unwrapped),30) + np.pi, 2*np.pi) - np.pi

# TODO flagging
###gaps longer than ROTATION are flagged
flag_10 = np.ceil(np.interp(ut_sci_10, ut[1:], np.diff(ut) > 80/3600.)).astype(np.uint8)
flag_15 = np.ceil(np.interp(ut_sci_15, ut[1:], np.diff(ut) > 80/3600.)).astype(np.uint8)
##self.synched_data[device]['FLAG'] = flag

#plt.figure()
#plt.plot(ut_sci_10, fixed_az_10, 'r.', label='fixed')
#plt.xlabel('UT')
#plt.savefig('fixedaz.png')

import pycfitsio as fits
fits.write(os.path.join(folder, 'fixaz.fits'), OrderedDict([
    ('10GHz', OrderedDict([('UT', ut_sci_10), ('AZ', fixed_az_10), ('FLAG', flag_10)])),
    ('15GHz', OrderedDict([('UT', ut_sci_15), ('AZ', fixed_az_15), ('FLAG', flag_15)]))
    ]))
