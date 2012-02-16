import exceptions
import pyfits
import pycfitsio as fits
import os
import numpy as np
from numpy.lib.recfunctions import append_fields
from collections import OrderedDict

REVCOUNTER_LABEL = {10:'REVCOUNTER_15GHZ', 15:'REVCOUNTER_10GHZ'}

WRAPPING_FIELDS_PI = ['HYBRIDHEADINGANGLE', 'HYBRIDYAWANGLE']
ROTATION = 80*2.e9
WRAPPING_FIELDS_360 = ['HEADING']

DEVICES = ['ANALOGCHANNELS', 'MAGNETOMETER', 'ASHTECH', 'GYRO_HID']

def find_clock_offsets_from_gpstime(cc, gpstime):
    """Computes computerClock offsets after restarts using gpstime
    
    Parameters
    ----------
    cc : ndarray
        computerClock array
    gpstime : ndarray
        gpstime array

    Returns:
    offsets : ndarray
        offsets to apply to computerclock to compensate for restarts
    """

    print('Find computerClock offsets')
    offsets_indices, = np.where(np.diff(cc)<0)
    print('Found clock jumps at indices %s, at relative position %s' % (str(offsets_indices), str(offsets_indices.astype(np.double)/len(cc)) ))
    offsets = [] 
    for j in offsets_indices:
        offsets.append(cc[j] + 2e9 * (gpstime[j+1] - gpstime[j]) - cc[j+1])
    return offsets

def apply_cc_offsets(cc, offsets):
    jumps, = np.where(np.diff(cc)<0)
    # apply offsets estimated with gpstime to computerclock of the revcounter
    if len(jumps) < len(offsets):
        print('Missing data in device')
    for index, offset in zip(jumps, offsets[:len(jumps)):
        cc[index+1:] += offset
    return cc

def read_data(base_folder='/home/zonca/COFE/data/sync_data', day='all', freq=10):
    gyro = fits.read(os.path.join(base_folder, 'servo', '%s.fits' % day), 'GYRO_HID')
    revcounter = fits.read(os.path.join(base_folder, 'servo', '%s.fits' % day), REVCOUNTER_LABEL[freq])
    data_rev = fits.read(os.path.join(base_folder, '%d' % freq, '%s.fits' % day), 0)['REV']
    return gyro, revcounter, data_rev

def create_science_computerclock(gyro, revcounter, data_rev):
    servo_range = remove_reset(revcounter['VALUE'], offsetsci=data_rev[0])

    # apply offsets to revcounter cc
    revcounter['computerClock'] = apply_cc_offsets(revcounter['computerClock'], offsets)

    # oversample revcounter cc and value to 140 Hz in order to interpolate over gaps
    uniform_rev_cc = np.arange(revcounter['computerClock'][servo_range][0], revcounter['computerClock'][servo_range][-1], (1/140.)*2e9, dtype=np.double)
    uniform_rev = np.interp( uniform_rev_cc, revcounter['computerClock'][servo_range].astype(np.double), revcounter['VALUE'][servo_range].astype(np.double))

    # create science data computer clock
    sci_cc = np.around(np.interp(data_rev, uniform_rev, uniform_rev_cc).astype(np.long)
    return sci_cc, offsets

def create_ut(gyro):
    """Create UT array from gpstime"""
    # check that gyro computerclock is already fixed
    assert np.all(np.diff(gyro['computerClock']) >= 0)
    # Fix gpstime to create the UT column
    ut = np.mod((gyro['GPSTIME'] + 15.)/3600., 24.)
    # remove single sample jumps
    good_ut = np.ones(len(ut),dtype=np.bool)
    good_ut[index_of_single_sample_jumps(ut)] = False

    # unwrap ut at 24 hours
    day_change_index, = np.where(np.diff(ut)<-23)
    assert len(day_change_index) == 1
    ut[day_change_index[0]+1:] += 24

    # get just the good samples
    utcc = gyro['computerClock'][good_ut]
    ut = ut[good_ut]
    return utcc, ut

def create_utscience(sci_file, offsets, utcc, ut, freq=10)
    """Create file with science data with fixed CC and UT
    see create_utservo
    """

    data = fits.read(sci_file)

    splitted_data['TIME'] = OrderedDict()
    for ch_n in range(16):
        ch_name = 'CHANNEL_%02d' % ch_n
        splitted_data[ch_name] = OrderedDict()
        for comp in 'TQU':
            splitted_data[ch_name][comp] = data[ch_name + comp]

def create_utservo(servo_file, offsets, utcc, ut)
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

    utservo = OrderedDict()
    for device in DEVICES:
        print('Processing ' + device)
        utservo[device] = fits.read(servo_file, device)
        utservo[device]['computerClock'] = apply_cc_offsets(utservo[device]['computerClock'], offsets)
        utservo[device]['UT'] = np.interp( utservo[device]['computerClock'], utcc, ut)

    filename = 'utservo.fits'
    print('Writing ' + filename)
    fits.write(filename, utservo)

def process_level1():
    gyro, revcounter, data_rev = read_data(base_folder='/home/zonca/COFE/data/sync_data', day='all', freq=10):
    offsets = find_clock_offsets_from_gpstime(gyro['computerClock'], gyro['GPSTIME'])
    utcc, ut = create_ut(gyro)
    create_utservo(servo_file, offsets, utcc, ut)
    create_utscience(sci_file, offsets, utcc, ut, freq=10)

class ServoSciSync(object):
    """Synchronizes servo and science data in a single fits file per day"""

    def __init__(self, base_folder = '/home/zonca/COFE/data/sync_data', day = '20110224', freq = 15, version=None):
        self.base_folder = base_folder
        self.day = day
        self.freq = freq
        self.synched_data = OrderedDict()
        self.version = version
        try:
            os.mkdir(os.path.join(self.base_folder, 'Level1'))
        except:
            pass
        self.out_folder = os.path.join(self.base_folder, 'Level1', '%s' % self.version)
        try:
            os.mkdir(self.out_folder)
        except:
            pass

    def load_data(self):
        self.servo = pyfits.open(os.path.join(self.base_folder, 'servo', '%s.fits' % self.day))
        self.devices = [ext.name for ext in self.servo[1:] if not ext.name.startswith('REV') and not ext.name.startswith('DEVICE') and not ext.name.startswith('TELESCOPE')]
                    
        
        self.data = pycfitsio.read(
                os.path.join(self.base_folder, 
                            '%d' % self.freq, 
                            '%s.fits' % self.day), 0)

        self.splitted_data = OrderedDict()
        self.splitted_data['TIME'] = OrderedDict()
        for ch_n in range(16):
            ch_name = 'CHANNEL_%02d' % ch_n
            self.splitted_data[ch_name] = OrderedDict()
            for comp in 'TQU':
                self.splitted_data[ch_name][comp] = self.data[ch_name + comp]

    def fix_counters(self):
        print('Fixing counters')

        self.counters = {}
        servo_count = self.servo[REVCOUNTER_LABEL[self.freq]].data.field(2)
        sci_count = self.data['REV']
        servo_range = remove_reset(servo_count, offsetsci=sci_count[0])
        #sci_range = remove_reset(sci_count)
        self.counters = {
                    'servo_range' : servo_range,
                    'servo' : servo_count[servo_range],
                    'servocc': self.servo[REVCOUNTER_LABEL[self.freq]].data['computerClock'][servo_range],
                    'sci_range' : None,
                    'sci' : sci_count #fix_counter_jumps(sci_count)
                 } 

    def find_clock_offsets_from_gpstime(self):
        print('Find computerClock offsets')
        device = 'GYRO_HID'
        cc = self.servo[device].data.field('computerClock')
        gpstime = self.servo[device].data.field('GPSTIME')
        self.offsets_indices, = np.where(np.diff(cc)<0)
        print('Found clock jumps at indices %s, at relative position %s' % (str(self.offsets_indices), str(self.offsets_indices.astype(np.double)/len(cc)) ))
        self.offsets = [] 
        for j in self.offsets_indices:
            self.offsets.append(cc[j] + 2e9 * (gpstime[j+1] - gpstime[j]) - cc[j+1])

        self.offsets_cc = cc[self.offsets_indices]  
        self.offsets_cc[1:] += np.cumsum(self.offsets[:-1])
            
    def sync_clock(self):
        print('SCI computer clock')
        assert np.all(np.diff(self.counters['servo']) >= 0)
        #make_monotonic(self.servo[REVCOUNTER_LABEL[self.freq]].data.field('computerClock'))
        cc = self.servo[REVCOUNTER_LABEL[self.freq]].data.field('computerClock')
        jumps, = np.where(np.diff(cc)<0)
        # apply offsets estimated with gpstime to computerclock of the revcounter
        assert len(jumps) == len(self.offsets)
        for index, offset in zip(jumps, self.offsets):
            cc[index+1:] += offset

        self.counters['fservocc'] = np.arange(cc[self.counters['servo_range']][0], cc[self.counters['servo_range']][-1], (1/140.)*2e9 )
        self.counters['fservo'] = np.interp( self.counters['fservocc'], cc[self.counters['servo_range']], self.counters['servo'])
        self.counters['flag'] = np.ceil(np.interp(self.counters['fservocc'], cc[self.counters['servo_range']][1:], np.diff(cc[self.counters['servo_range']])>ROTATION))
        self.synched_data['TIME'] = OrderedDict()
        self.synched_data['TIME']['computerClock'] = np.around(np.interp(self.counters['sci'], self.counters['fservo'], self.counters['fservocc'])).astype(np.int64)
        #self.synched_data['computerClock'] = np.around(np.interp(self.counters['sci'], self.counters['servo'], 
        #            cc[self.counters['servo_range']])).astype(np.int64)
        self.splitted_data['TIME']['computerClock'] = self.synched_data['TIME']['computerClock']
        self.synched_data['TIME']['norevcountflag'] = np.ceil(np.interp(self.synched_data['TIME']['computerClock'], self.counters['fservocc'], self.counters['flag'])).astype(np.uint8)
        self.splitted_data['TIME']['norevcountflag'] = self.synched_data['TIME']['norevcountflag'] 

    def fix_devices_cc(self, cc_int):
        cc = cc_int.astype(np.double)
        jumps, = np.where(np.diff(cc)<0)
        # apply offsets estimated with gpstime to computerclock of each device
        offsets = self.offsets
        if len(jumps) < len(self.offsets):
            print('Missing data in device')
            offsets = self.offsets[:len(jumps)]
        for index, offset in zip(jumps, offsets):
            cc[index+1:] += offset
        return cc

    def sync_devices(self):
        print('Synching devices')
        for device in self.devices:
            print(device)
            ext = self.servo[device]
            self.synched_data[device] = OrderedDict()
            cc = self.fix_devices_cc(ext.data.field('computerClock'))

            #gaps longer than ROTATION are flagged
            flag = np.ceil(np.interp(self.splitted_data['TIME']['computerClock'],cc[1:], np.diff(cc) > ROTATION)).astype(np.uint8)
            self.synched_data[device]['FLAG'] = flag

            assert np.all(np.diff(cc) >= 0)
            for col in ext.columns[1:]:
                if col.name in WRAPPING_FIELDS_PI:
                    print('UNWRAP ' + col.name)
                    valid = (ext.data[col.name] < np.pi) & (ext.data[col.name] > -np.pi)

                    #unwrap the heading angle 
                    wraps = np.diff(col.array[valid]) < - .95 * 2 * np.pi #5% tolerance
                    unwrapped = col.array[valid].copy()
                    unwrapped[1:] += np.cumsum(wraps) * np.pi * 2

                    #fix time gaps
                    typical_revlength = np.median(np.diff(cc[wraps]))
                    h_jumps = np.diff(cc[valid]) > typical_revlength * .8
                    h_jumps_scaled = h_jumps.astype(np.double) 
                    h_jumps_scaled[h_jumps] *= np.round(np.diff(cc[valid])[h_jumps]/typical_revlength)
                    unwrapped[1:] += np.cumsum(h_jumps_scaled) * np.pi * 2 

                    #interpolate and reset to -pi pi
                    self.synched_data[device][col.name] = np.mod(np.interp(self.splitted_data['TIME']['computerClock'], cc[valid], unwrapped) + np.pi, 2*np.pi) - np.pi
                else:
                    try:
                        self.synched_data[device][col.name] = np.interp(self.splitted_data['TIME']['computerClock'], cc, 
                                    ext.data.field(col.name))
                    except exceptions.ValueError:
                        print('SKIPPING %s, no samples in range' % '_'.join([ext.name, col.name]))
                        self.synched_data['TIME']['REVCHECK'] = np.interp(self.splitted_data['TIME']['computerClock'],
                            self.servo[REVCOUNTER_LABEL[self.freq]].data.field('computerClock')[self.counters['servo_range']],
                            self.servo[REVCOUNTER_LABEL[self.freq]].data.field('value')[self.counters['servo_range']]
                            )
    #self.synched_data['TIME']['UT'] = cctout(self.splitted_data['TIME']['computerClock'], self.synched_data['GYRO_HID']['GPSTIME'])
            if device == 'GYRO_HID':
                #Fix gpstime to create the UT column
                ut = np.mod((self.servo[device].data['GPSTIME'] + 15.)/3600., 24.)
                #remove single sample jumps
                ut_mask = np.zeros(len(ut),dtype=np.bool)
                ut_mask[index_of_single_sample_jumps(ut)] = True
                ut = np.ma.MaskedArray(ut, mask=ut_mask)
                day_change_index, = np.where(np.diff(ut)<-23)
                assert len(day_change_index) == 1
                ut[day_change_index[0]+1:] += 24

                # Leave ut in original sampling rate
                self.utcc = np.ma.MaskedArray(cc, mask=ut.mask).compressed()
                self.ut = ut.compressed()

                #self.synched_data['TIME']['UT'] = cctout(self.splitted_data['TIME']['computerClock'], self.synched_data['GYRO_HID']['GPSTIME'])
                self.synched_data['TIME']['UT'] = np.interp(self.splitted_data['TIME']['computerClock'],  self.utcc, self.ut)
                self.splitted_data['TIME']['UT'] = self.synched_data['TIME']['UT']
                self.offsets_ut = np.interp(self.offsets_cc,  self.utcc, self.ut)
                np.savetxt(self.out_folder + '/cc_offsets_%dGHz_cc-ut-off.txt' % self.freq, (self.offsets_cc, self.offsets_ut, self.offsets))

    def write(self):
        filename = os.path.join(self.out_folder,'%s_%dGHz_v%s' % (self.day, self.freq, self.version))
        print('Writing %s' % filename)
        #np.save(filename.replace('.fits','_data.npy', self.data))
        #np.save(filename.replace('.fits','_servo.npy', self.synched_data))
        pycfitsio.write(filename +  '_data.fits', self.splitted_data)
        pycfitsio.write(filename + '_servo.fits', self.synched_data)

    def write_devices_with_ut(self):
        """Devices have already a CC column fixed, so we just interpolate the UT to that sampling"""

        #free memory
        try:
            del self.splitted_data
        except:
            pass

        filename = os.path.join(self.out_folder,'utservo_%s_v%s.fits' % (self.day, self.version))
        self.utservo = OrderedDict()
        for device in self.devices:
            self.utservo[device] = OrderedDict()
            self.utservo[device]['UT'] = np.interp(self.fix_devices_cc(self.servo[device].data['computerClock']),  self.utcc, self.ut)
            for colname in self.servo[device].data.dtype.names:
                self.utservo[device][colname] = np.array(self.servo[device].data[colname])

        print('Writing ' + filename)
        pycfitsio.write(filename, self.utservo)

    def process(self):
        self.load_data()
        self.fix_counters()
        self.find_clock_offsets_from_gpstime()
        self.sync_clock()
        self.sync_devices()

    def run(self):
        self.process()
        self.write()
        if self.freq == 10: # needed just once
            self.write_devices_with_ut()
