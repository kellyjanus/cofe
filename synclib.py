import exceptions
import pyfits
import pycfitsio
import os
import numpy as np
from numpy.lib.recfunctions import append_fields
from collections import OrderedDict

REVCOUNTER_LABEL = {10:'REVCOUNTER_15GHZ', 15:'REVCOUNTER_10GHZ'}

WRAPPING_FIELDS_PI = ['HYBRIDHEADINGANGLE', 'HYBRIDYAWANGLE']
ROTATION = 80*2.e9
WRAPPING_FIELDS_360 = ['HEADING']

def cctotime(c):
    return (c-c[0])*2.e-9

def cctout(cc, gpstime):
    """cc and gpstime must be already synched"""
    utc = np.mod(((gpstime+15.)/3600.), 24)
    ut = utc[0]+(((cc - cc[0])*2e-9)/3600.)
    return ut

def fix_counter_jumps(d):
    """Removes jumps by linear fitting and removing points further than a predefined threshold, then linearly interpolates to the full array"""
    THRESHOLD = 3
    t = np.arange(len(d))
    ar, br = np.polyfit(t,d,1)
    lin_d = np.polyval([ar, br], t)
    valid = np.abs(d - lin_d) < THRESHOLD
    d_fix = np.interp(t, t[valid], d[valid], left=d[valid][0]-1, right=d[valid][-1]+1)
    return d_fix.astype(np.int)

def fix_counter_jumps_diff(d):
    """Removes 1 sample jumps by checking the sample to sample difference"""
    THRESHOLD = 5
    t = np.arange(len(d))
    fix = np.abs(np.diff(d))<THRESHOLD
    lin_d = d[fix]
    d_fix = np.interp(t, t[fix], lin_d, left=lin_d[0]-1, right=lin_d[-1]+1)
    return d_fix.astype(np.int)

def remove_reset(d, offsetsci=None):
    """Gets longest time period between resets"""
    # first check for 20bit jumps
    jump20bit_indices, = np.where(np.logical_and(np.diff(d) < -2**20*.9, np.diff(d) > -2**20*1.1))
    print('20bit jumps at:')
    print(jump20bit_indices)
    for j in jump20bit_indices:
        d[j+1:] += 2**20
    if not offsetsci is None:
        d += 2**20 * np.round((offsetsci - d[0])/2**20)
    print('remove dip')
    start,=np.where(np.diff(d) < -500000)
    stop,=np.where(np.diff(d) > 500000)
    if len(start)>0:
        d[start[1]-1:stop[1]+1] = -1
    reset_indices, = np.where(np.diff(d) < -300000)
    real_reset = []
    for i in reset_indices:
        if not (1.5e6<i<1.6e6): 
        #if abs(d[i+2]-d[i]) >= abs(d[i+1]-d[i]) and not (1.5e6<i<1.6e6): 
            #it is a single point, not a reset
            real_reset.append(i)
    reset_indices = np.array(real_reset)

    print('reset jumps at:')
    print(reset_indices)
    sections_boundaries = np.concatenate([[0], reset_indices +1 ,[len(d)]])
    sections_lengths = np.diff(sections_boundaries)
    max_len = sections_lengths.argmax()
    s = np.zeros(len(d), dtype=np.bool)
    start = np.max([150000, sections_boundaries[max_len]])
    stop = np.min([sections_boundaries[max_len+1], d.searchsorted(4865400)])
    s[start+1:stop-1] = True
    s[d==-1] = False
    s[d==0] = False
    single_sample_jump, = np.where(np.diff(d[:stop-1])<0)
    assert np.all(np.diff(single_sample_jump)>2)
    s[single_sample_jump + 1] = False
    print('selecting section with no resets from index %d to index %d %.2f %% of the total length' % (start, stop, (stop-start)*100./len(d)))
    return s

def index_of_single_sample_jumps(d):
    jumps,=np.where(np.abs(np.diff(d))>.005 )
    single_jumps = jumps[np.diff(jumps)==1] + 1
    return single_jumps

def make_monotonic(d):
    jumps, = np.where(np.diff(d)<0)
    print('Fixing jumps at indices %s, at relative position %s' % (str(jumps), str(jumps.astype(np.double)/len(d)) ))
    for j in jumps:
        d[j+1:] += d[j] - d[j+1]

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
                    'sci_range' : None,
                    'sci' : sci_count #fix_counter_jumps(sci_count)
                 } 

    def find_clock_offsets_from_gpstime(self):
        print('Find computerClock offsets')
        device = 'GYRO_HID'
        cc = self.servo[device].data.field('computerClock')
        gpstime = self.servo[device].data.field('GPSTIME')
        jumps, = np.where(np.diff(cc)<0)
        print('Found clock jumps at indices %s, at relative position %s' % (str(jumps), str(jumps.astype(np.double)/len(cc)) ))
        self.offsets = [] 
        for j in jumps:
            self.offsets.append(cc[j] + 2e9 * (gpstime[j+1] - gpstime[j]) - cc[j+1])
            
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

        self.counters['fservocc'] = np.arange( cc[self.counters['servo_range']][0], cc[self.counters['servo_range']][-1], (1/40.)*2e9 )
        self.counters['fservo'] = np.interp( self.counters['fservocc'], cc[self.counters['servo_range']], self.counters['servo'])
        self.synched_data['TIME'] = OrderedDict()
        self.synched_data['TIME']['computerClock'] = np.around(np.interp(self.counters['sci'], self.counters['fservo'], self.counters['fservocc'])).astype(np.int64)
        #self.synched_data['computerClock'] = np.around(np.interp(self.counters['sci'], self.counters['servo'], 
        #            cc[self.counters['servo_range']])).astype(np.int64)
        self.splitted_data['TIME']['computerClock'] = self.synched_data['TIME']['computerClock']

    def sync_devices(self):
        print('Synching devices')
        for device in self.devices:
            print(device)
            ext = self.servo[device]
            self.synched_data[device] = OrderedDict()
            cc = ext.data.field('computerClock')

            jumps, = np.where(np.diff(cc)<0)
        # apply offsets estimated with gpstime to computerclock of each device
            if len(jumps) > 0:
                offsets = self.offsets
                if len(jumps) < len(self.offsets):
                    print('Missing data in %s' % device)
                    offsets = self.offsets[:len(jumps)]
                for index, offset in zip(jumps, offsets):
                    cc[index+1:] += offset

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


                    #self.synched_data['TIME']['UT'] = cctout(self.splitted_data['TIME']['computerClock'], self.synched_data['GYRO_HID']['GPSTIME'])
                    self.synched_data['TIME']['UT'] = np.interp(self.splitted_data['TIME']['computerClock'], np.ma.MaskedArray(cc, mask=ut.mask).compressed(), ut.compressed()) 
                    self.splitted_data['TIME']['UT'] = self.synched_data['TIME']['UT']

    def write(self):
        folder = os.path.join(self.base_folder, 'Level1', '%s' % self.version)
        try:
            os.mkdir(folder)
        except:
            pass
        filename = os.path.join(folder,'%s_%dGHz_v%s' % (self.day, self.freq, self.version))
        print('Writing %s' % filename)
        #np.save(filename.replace('.fits','_data.npy', self.data))
        #np.save(filename.replace('.fits','_servo.npy', self.synched_data))
        pycfitsio.write(filename +  '_data.fits', self.splitted_data)
        pycfitsio.write(filename + '_servo.fits', self.synched_data)


    def run(self):
        self.load_data()
        self.fix_counters()
        self.find_clock_offsets_from_gpstime()
        self.sync_clock()
        self.sync_devices()
        self.write()
