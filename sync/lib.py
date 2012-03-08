import exceptions
import pycfitsio as fits
import os
import numpy as np
from collections import OrderedDict

from utils import *

REVCOUNTER_LABEL = {10:'REVCOUNTER_15GHZ', 15:'REVCOUNTER_10GHZ'}

WRAPPING_FIELDS_PI = ['HYBRIDHEADINGANGLE', 'HYBRIDYAWANGLE']
ROTATION = 80/2e-9
WRAPPING_FIELDS_360 = ['HEADING']

DEVICES = ['ANALOGCHANNELS', 'MAGNETOMETER', 'ASHTECH', 'GYRO_HID']

def find_clock_offsets_from_gpstime(cc, gpstime):
    """Computes computerClock offsets after restarts using gpstime

    all negative jumps in computerClock mean a restart.
    if the restart is between sample j and j+1, we compute
    the needed offset as:
    cc[j] - cc[j+1] + (gpstime[j+1] - gpstime[j])/2e-9
    where:
    cc[j+1] - cc[j] is the current jump in computerClock
    (gpstime[j+1] - gpstime[j]) /2e-9 is the correct jump computed from gpstime    so we actually remove the current jump and add back
    the correct jump converted from gpstime
    
    Parameters
    ----------
    cc : ndarray
        computerClock array
    gpstime : ndarray
        gpstime array

    Returns
    -------
    offsets : ndarray
        offsets to apply to computerclock to compensate for restarts
    """

    print('Find computerClock offsets')
    offsets_indices, = np.where(np.diff(cc)<0)
    print('Found clock jumps at indices %s, at relative position %s' % (str(offsets_indices), str(offsets_indices.astype(np.double)/len(cc)) ))
    offsets = [] 
    for j in offsets_indices:
        offsets.append(cc[j] + (gpstime[j+1] - gpstime[j])/2e-9 - cc[j+1])
    cc = apply_cc_offsets(cc, offsets)
    return offsets

def apply_cc_offsets(cc, offsets):
    """Blindly apply offsets computed from gpstime to cc

    offsets are blindly summed to the negative jumps
    in the cc array.

    Parameters
    ----------
    cc : ndarray
        computerClock to be corrected for restarts
    offsets : ndarray
        precomputed offsets

    Returns
    -------
    cc : ndarray
        corrected computerClock
    """
    jumps, = np.where(np.diff(cc)<0)
    # apply offsets estimated with gpstime to computerclock of the revcounter
    if len(jumps) < len(offsets):
        print('Missing data in device')
    for index, offset in zip(jumps, offsets[:len(jumps)]):
        cc[index+1:] += offset
    return cc

def create_science_computerclock(gyro, revcounter, data_rev, offsets):
    """Syncronizes science and gyro using the revcounter

    the revcounter has several gaps, so we first interpolate it uniformly
    at 140 Hz and then we interpolate the uniform revcounter computerclock
    to the science sampling rate using the revcounters

    Parameters
    ----------
    gyro : OrderedDict
        gyro data
    revcounter : OrderedDict
        revcounter data
    data_rev : ndarray
        revcounter of the scientific channel
    offsets : ndarray
        cc offsets from gpstime

    Returns
    -------
    sci_cc : ndarray
        syncronized scientific computerclock
    """
    # servo revcounter is cleaned up from restarts and jumps
    servo_range = remove_reset(revcounter['VALUE'], offsetsci=data_rev[0])

    # apply offsets to revcounter cc
    revcounter['COMPUTERCLOCK'] = apply_cc_offsets(revcounter['COMPUTERCLOCK'], offsets)

    # oversample revcounter cc and value to 140 Hz in order to interpolate over gaps
    uniform_rev_cc = np.arange(revcounter['COMPUTERCLOCK'][servo_range][0], revcounter['COMPUTERCLOCK'][servo_range][-1], (1/1000.)/2.e-9, dtype=np.double)
    uniform_rev = np.interp( uniform_rev_cc, revcounter['COMPUTERCLOCK'][servo_range].astype(np.double), revcounter['VALUE'][servo_range].astype(np.double))

    flag = np.ceil(np.interp(uniform_rev_cc, revcounter['COMPUTERCLOCK'][1:], np.diff(revcounter['COMPUTERCLOCK'])>ROTATION))


    # create science data computer clock
    sci_cc = np.around(np.interp(data_rev.astype(np.double), uniform_rev.astype(np.double), uniform_rev_cc.astype(np.double)).astype(np.long))

    norevcountflag = np.ceil(np.interp(sci_cc, uniform_rev_cc, flag)).astype(np.uint8)

    return sci_cc, norevcountflag

def create_ut_from_gpstime(gpstime):
    # Fix gpstime to create the UT column
    fixed = np.mod((gpstime + 15.)/3600., 24.)
    # remove single sample jumps
    good = np.ones(len(fixed),dtype=np.bool)
    good[index_of_single_sample_jumps(fixed)] = False

    # get just the good samples
    gpstime_index = np.arange(len(gpstime))[good]

    # interpolate back to original length
    fixed = np.interp(np.arange(len(gpstime)), gpstime_index, fixed[good])

    # unwrap ut at 24 hours
    day_change_index, = np.where(np.diff(fixed)<-23)
    assert len(day_change_index) == 1
    fixed[day_change_index[0]+1:] += 24
    return fixed

def create_ut_from_cc(gyro):
    """Create UT array from computerclock

    UT is defined as UT hour of the first day,
    e.g. 10 is 10am of the first day
    and monotonically increasing after 24,
    so 27.5 is 3:30am of the second day.
    
    Parameters
    ----------
    gyro : OrderedDict
        gyro data
        
    Returns
    -------
    utcc : ndarray
        computerclock of ut array
    ut : ndarray
        ut array
    """
    # check that gyro computerclock is already fixed
    assert np.all(np.diff(gyro['COMPUTERCLOCK']) >= 0)
    ut = np.mod((gyro['GPSTIME'] + 15.)/3600., 24.)

    utcc = gyro['COMPUTERCLOCK']
    ut = ut[0] + (gyro['COMPUTERCLOCK'] - gyro['COMPUTERCLOCK'][0]) * 2e-9 / 3600.

    return utcc, ut

def create_utscience(sci_file, gyro, revcounter, offsets, utcc, ut, freq):
    """Create file with science data with fixed CC and UT
    see create_utservo
    """

    data = fits.read(sci_file)

    splitted_data = OrderedDict()
    splitted_data['TIME'] = OrderedDict()
    for ch_n in range(16):
        ch_name = 'CH%d_' % ch_n
        splitted_data[ch_name] = OrderedDict()
        for comp in 'TQU':
            splitted_data[ch_name][comp] = data[ch_name + comp]
    splitted_data['TIME']['COMPUTERCLOCK'], splitted_data['TIME']['NOREVCOUNTFLAG']= create_science_computerclock(gyro, revcounter, data['REV'], offsets)
    splitted_data['TIME']['UT'] = np.interp(splitted_data['TIME']['COMPUTERCLOCK'], utcc, ut)

    filename = '%s_%dGHz_data.fits' % (os.path.basename(sci_file).split('.')[0], freq)
    print('Writing ' + filename)
    fits.write(filename, splitted_data)

    return splitted_data

def create_utservo(servo_file, offsets, utcc, ut):
    """Create file with servo data with fixed CC and UT

    Parameters
    ----------
    servo_file : str
        filename of servo data
    offsets : ndarray
        CC offsets computed from gpstime, see find_clock_offsets_from_gpstime
    utcc : ndarray
        CC array related to UT
    ut : ndarray
        UT array

    Returns
    -------
    writes utservo.fits file to disk
    """
    print("Create UT timestamped servo data")

    utservo = OrderedDict()
    for device in DEVICES:
        print('Processing ' + device)
        utservo[device] = fits.read(servo_file, device)
        utservo[device]['COMPUTERCLOCK'] = apply_cc_offsets(utservo[device]['COMPUTERCLOCK'], offsets)
        utservo[device]['UT'] = np.interp( utservo[device]['COMPUTERCLOCK'], utcc, ut)

    filename = 'utservo.fits'
    print('Writing ' + filename)
    fits.write(filename, utservo)
    return utservo

def create_sync_servo(servo_file, offsets, utcc, ut, sci_cc, freq):
    """Create synchronized servo data

    Parameters
    ----------
    servo_file : srt
        filename of servo data
    offsets : ndarray
        CC offsets computed from gpstime, see find_clock_offsets_from_gpstime
    utcc : ndarray
        CC array related to UT
    ut : ndarray
        UT array
    sci_cc : ndarray
        CC array of scientific data
    freq : int
        frequency
    """

    print("Create synchronized servo data to %dGHz" % freq)
    filename = '%s_%dGHz_servo.fits' % (os.path.basename(servo_file).split('.')[0], freq)
    print('Writing ' + filename)
    f = fits.create(filename)
    ext = OrderedDict()
    ext['COMPUTERCLOCK'] = sci_cc
    ext['UT'] = np.interp(sci_cc, utcc, ut)
    f.write_HDU('TIME', ext)
    for device in DEVICES:
        print('Processing ' + device)
        raw_data = fits.read(servo_file, device)
        ext = OrderedDict()
        cc = apply_cc_offsets(raw_data['COMPUTERCLOCK'], offsets)
        for colname, colarray in raw_data.iteritems():
            if colname != 'COMPUTERCLOCK':
                print('Column ' + colname)
                ext[colname] = np.interp(sci_cc, cc, colarray)
        f.write_HDU(device, ext)
    f.close()
    
def fix_gyro_cc(gyro, ut):
    return gyro['COMPUTERCLOCK'][0] + (ut-ut[0])*3600/2e-9
    
def process_level1(base_folder='/COFE', day='all', use_cc=True):
    """Full processing to produce Level1 data
    
    Parameters
    ----------
    base_folder : str
        path to data
    day : str
        day to be processed
    freq : int
        frequency
    """
    gyro = fits.read(os.path.join(base_folder, 'servo', '%s.fits' % day), 'GYRO_HID')
    offsets = find_clock_offsets_from_gpstime(gyro['COMPUTERCLOCK'], gyro['GPSTIME'])
    if use_cc:
        utcc, ut = create_ut_from_cc(gyro)
    else:
        ut = create_ut_from_gpstime(gyro['GPSTIME'])
        gyro['COMPUTERCLOCK'] = fix_gyro_cc(gyro, ut)
        utcc = gyro['COMPUTERCLOCK']
    servo_file = os.path.join(base_folder,'servo','%s.fits' % day)
    create_utservo(servo_file, offsets, utcc, ut)
    for freq in [10, 15]:
        revcounter = fits.read(os.path.join(base_folder, 'servo', '%s.fits' % day), REVCOUNTER_LABEL[freq])
        sci = create_utscience(os.path.join(base_folder,str(freq),'%s.fits'%day), gyro, revcounter, offsets, utcc, ut, freq)
        create_sync_servo(servo_file, offsets, utcc, ut, sci['TIME']['COMPUTERCLOCK'], freq)

if __name__ == '__main__':
    process_level1()
    #base_folder='/COFE'; day='all'
    #gyro = fits.read(os.path.join(base_folder, 'servo', '%s.fits' % day), 'GYRO_HID')
    #offsets = find_clock_offsets_from_gpstime(gyro['COMPUTERCLOCK'], gyro['GPSTIME'])
    #use_cc = True

    #servo_file = os.path.join(base_folder,'servo','%s.fits' % day)
    #utservo = create_utservo(servo_file, offsets, utcc, ut)
    #freq = 10
    #revcounter = fits.read(os.path.join(base_folder, 'servo', '%s.fits' % day), REVCOUNTER_LABEL[freq])
    #sci = create_utscience(os.path.join(base_folder,str(freq),'%s.fits'%day), gyro, revcounter, offsets, utcc, ut, freq)
    ##create_sync_servo(servo_file, offsets, utcc, ut, sci_cc, freq)
    #ra=slice(866653+1500,869045-700+1)
    #figure()
    #plot(sci['TIME']['UT'][ra],sci['CH1_']['T'][ra])
    #from smooth import smooth
    #figure()
    #plot(np.interp( sci['TIME']['UT'][ra], utservo['GYRO_HID']['UT'], utservo['GYRO_HID']['HYBRIDHEADINGANGLE'] ), sci['CH1_']['T'][ra])
    #plot(smooth(np.interp( sci['TIME']['UT'][ra], utservo['GYRO_HID']['UT'], utservo['GYRO_HID']['HYBRIDHEADINGANGLE'] ),30), sci['CH1_']['T'][ra])
