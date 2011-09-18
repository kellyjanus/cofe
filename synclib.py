import pyfits
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
    reset_indices, = np.where(np.diff(d) < -1000)
    sections_boundaries = np.concatenate([[0], reset_indices +1 ,[len(d)]])
    sections_lengths = np.diff(sections_boundaries)
    max_len = sections_lengths.argmax()
    return slice(sections_boundaries[max_len],sections_boundaries[max_len+1])

class ServoSciSync(object):
    """Synchronizes servo and science data in a single fits file per day"""

    def __init__(self, base_folder = '/home/zonca/COFE/data/sync_data', day = '20110224', freq = [10, 15]):
        self.base_folder = base_folder
        self.day = day
        self.freq = freq
        try:
            os.mkdir(os.path.join(self.base_folder, 'Level1'))
        except:
            pass

    def load_data(self):
        self.servo = pyfits.open(os.path.join(self.base_folder, 'Servo_Test_Data', '%s.fits' % self.day))
        self.devices = [ext.name for ext in self.servo[1:] if not ext.name.startswith('REV')]
        
        self.data = {}
        for freq in self.freq:
            pyf = pyfits.open(
                    os.path.join(self.base_folder, 
                                '%d_GHz_Data' % freq, 
                                '%s.fits' % self.day))
            pyf.readall()
            self.data[freq] = pyf[1].columns.data

    def fix_counters(self):

        self.counters = {}
        for freq in self.freq:
            servo_count = self.servo['REVCOUNTER_%dGHZ' % freq].data.field(2)
            servo_range = remove_reset(servo_count)
            sci_count = self.data[freq][0].array
            #sci_range = remove_reset(sci_count)
            self.counters[freq] = {
                        'servo_range' : servo_range,
                        'servo' : fix_counter_jumps_diff(servo_count[servo_range]),
                        'sci_range' : None,
                        'sci' : fix_counter_jumps(sci_count)
                     } 

    def sync_clock(self):
        self.clock = {}
        for freq in self.freq:
            self.data[freq][0].array[:] = self.counters[freq]['sci']
            self.clock[freq] =np.around(np.interp(self.counters[freq]['sci'], self.counters[freq]['servo'], self.servo['REVCOUNTER_%dGHZ' % freq].data.field('computerClock')[self.counters[freq]['servo_range']])).astype(np.int64)
            newcol = pyfits.Column(
                name = 'computerClock',
                format = 'K',
                array = self.clock[freq]
                )
            self.data[freq].append(newcol)

    def sync_devices(self):
        for freq in self.freq:
            for device in self.devices:
                ext = self.servo[device]
                for col in ext.columns[1:]:
                    newcol = pyfits.Column(
                                name = '_'.join([ext.name, col.name]),
                                format = col.format,
                                array = np.interp(self.clock[freq], ext.data.field('computerClock')[self.counters[freq]['servo_range']], ext.data.field(col.name)[self.counters[freq]['servo_range']])
                                )
                    self.data[freq].append(newcol)
            self.data[freq].append(
                    pyfits.Column(
                        name='REVCHECK',
                        format = 'K',
                        array = np.interp(self.clock[freq],
                            self.servo['REVCOUNTER_%dGHZ' % freq].data.field('computerClock')[self.counters[freq]['servo_range']],
                            self.servo['REVCOUNTER_%dGHZ' % freq].data.field('value')[self.counters[freq]['servo_range']]
                            )
                        )
                    )
    
    def write(self):
        for freq in self.freq:
            filename = os.path.join(self.base_folder, 'Level1','%s_%dGHz.fits' % (self.day, freq))
            print('Writing %s' % filename)
            table = pyfits.new_table(self.data[freq])
            table.name = 'DATA'
            table.writeto(filename,clobber=True)


    def run(self):
        self.load_data()
        self.fix_counters()
        self.sync_clock()
        self.sync_devices()
        self.write()

