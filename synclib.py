import pyfits
import pycfitsio
import os
import numpy as np

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

def remove_reset(d):
    """Gets longest time period between resets"""
    reset_indices, = np.where(np.diff(d) < -100000)
    sections_boundaries = np.concatenate([[0], reset_indices +1 ,[len(d)]])
    sections_lengths = np.diff(sections_boundaries)
    max_len = sections_lengths.argmax()
    return slice(sections_boundaries[max_len],sections_boundaries[max_len+1])

class ServoSciSync(object):
    """Synchronizes servo and science data in a single fits file per day"""

    def __init__(self, base_folder = '/home/zonca/COFE/data/sync_data', day = '20110224', freq = 15):
        self.base_folder = base_folder
        self.day = day
        self.freq = freq
        try:
            os.mkdir(os.path.join(self.base_folder, 'Level1'))
        except:
            pass

    def load_data(self):
        self.servo = pyfits.open(os.path.join(self.base_folder, 'servo', '%s.fits' % self.day))
        self.devices = [ext.name for ext in self.servo[1:] if not ext.name.startswith('REV') and not ext.name.startswith('DEVICE') and not ext.name.startswith('TELESCOPE')]
                    
        
        self.data = {}
        self.data, self.data_header = pycfitsio.read(
                os.path.join(self.base_folder, 
                            '%d' % self.freq, 
                            '%s.fits' % self.day), 0)

    def fix_counters(self):

        self.counters = {}
        servo_count = self.servo['REVCOUNTER_%dGHZ' % freq].data.field(2)
        servo_range = remove_reset(servo_count)
        sci_count = self.data[0].array
        #sci_range = remove_reset(sci_count)
        self.counters[freq] = {
                    'servo_range' : servo_range,
                    'servo' : fix_counter_jumps_diff(servo_count[servo_range]),
                    'sci_range' : None,
                    'sci' : fix_counter_jumps(sci_count)
                 } 

    def sync_clock(self):
        self.clock = {}
        self.data[0].array[:] = self.counters[freq]['sci']
        self.clock =np.around(np.interp(self.counters['sci'], self.counters['servo'], self.servo['REVCOUNTER_%dGHZ' % freq].data.field('computerClock')[self.counters['servo_range']])).astype(np.int64)
        newcol = pyfits.Column(
            name = 'computerClock',
            format = 'K',
            array = self.clock[freq]
            )
        self.data.append(newcol)

    def sync_devices(self):
            for device in self.devices:
                ext = self.servo[device]
                for col in ext.columns[1:]:
                    newcol = pyfits.Column(
                                name = '_'.join([ext.name, col.name]),
                                format = col.format,
                                array = np.interp(self.clock, ext.data.field('computerClock')[self.counters['servo_range']], ext.data.field(col.name)[self.counters['servo_range']])
                                )
                    self.data.append(newcol)
            self.data.append(
                    pyfits.Column(
                        name='REVCHECK',
                        format = 'K',
                        array = np.interp(self.clock[self.freq],
                            self.servo['REVCOUNTER_%dGHZ' % self.freq].data.field('computerClock')[self.counters[self.freq]['servo_range']],
                            self.servo['REVCOUNTER_%dGHZ' % self.freq].data.field('value')[self.counters[self.freq]['servo_range']]
                            )
                        )
                    )
    
    def write(self):
        filename = os.path.join(self.base_folder, 'Level1','%s_%dGHz.fits' % (self.day, self.freq))
        print('Writing %s' % filename)
        table = pyfits.new_table(self.data)
        table.name = 'DATA'
        table.writeto(filename,clobber=True)


    def run(self):
        self.load_data()
        self.fix_counters()
        self.sync_clock()
        self.sync_devices()
        self.write()

