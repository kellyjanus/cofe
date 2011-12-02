import matplotlib.pyplot as plt
import sys
sys.path.append('../code')

import synclib
from synclib import *

freq = 10
self = synclib.ServoSciSync(base_folder ='/COFE', day = 'all', freq = freq, version='debug')
self.load_data()
plt.figure()
plt.title('%d GHz REV' % freq)
plt.plot(self.data['REV'],'.',label='SCI REV before fix')
plt.plot(self.servo[REVCOUNTER_LABEL[self.freq]].data.field(2) ,'.',label='SERVO REV before fix')
self.fix_counters()
plt.plot(self.counters['sci'],'.',label='SCI REV after fix')
plt.plot(self.counters['servo'],'.',label='SERVO REV after fix')
plt.legend(loc=0); plt.grid(); plt.xlabel('Sample index')

plt.figure()
plt.title('Computer clocks before making monotonic')
for device in self.devices:
    plt.plot(self.servo[device].data.field('computerClock'), label = device)

self.find_clock_offsets_from_gpstime()
self.sync_clock()
plt.plot(self.splitted_data['TIME']['computerClock'], label = 'SCI')
plt.legend(loc=0); plt.grid()

self.sync_devices()
plt.figure()
plt.plot(self.synched_data['GYRO_HID']['GPSTIME'])
plt.grid()

plt.figure()
plt.plot(self.splitted_data['TIME']['computerClock'], label = 'SCI')
plt.title('Computer clocks after making monotonic')
for device in self.devices:
    plt.plot(self.servo[device].data.field('computerClock'), label = device)
plt.legend(loc=0); plt.grid()

plt.figure()
plt.title('%d GHz REV' % freq)
plt.plot(self.splitted_data['TIME']['computerClock'], self.counters['sci'],'r.',label='SCI REV')
plt.plot(self.synched_data['TIME']['computerClock'], self.synched_data['TIME']['REVCHECK'],'.', label='SERVO REV')
plt.legend(); plt.grid(); plt.xlabel('computerClock')

plt.figure()
plt.title('%d GHz MAG' % freq)
plt.plot(cctotime(self.splitted_data['TIME']['computerClock']), np.degrees(np.arctan2(self.splitted_data['CHANNEL_15']['T'], self.splitted_data['CHANNEL_14']['T'])),label='SCI')
plt.plot(cctotime(self.synched_data['TIME']['computerClock']), np.degrees(np.arctan2(self.synched_data['ANALOGCHANNELS']['CHANNEL_31'], self.synched_data['ANALOGCHANNELS']['CHANNEL_30'])),label='ANALOGCH')
plt.legend(); plt.grid(); plt.xlabel('hours')

#self.write()
##self.run()
plt.show()
