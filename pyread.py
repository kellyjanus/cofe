from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt

NCHAN = 16
SEC_PER_REV = 256
ENC_START_TRIGGER = 15

filename = '15475300.dat'
filename = '../01073800.dat'
filename = '../18402800.dat'

# the data type defines a new type of array, like an array of structures
channels_labels = ['ch%d' % i for i in range(NCHAN)]
cofedtype=np.dtype( [(ch,np.uint16) for ch in channels_labels] + [('enc',np.uint16)]+[('dummy',np.uint16)]  + [('rev%d' % i,np.uint16) for i in range(3)])

# the memory map works as an array but directly on the data on disk (if mode is not r, the file on disk is modified by array operations)
data_on_disk=np.memmap(filename,dtype=cofedtype,mode='r')

# IDL 16 bits are from 0 to 2**16 instead python is -2**15 + 2**15
# data are accessible as d['ch0'] ... d['ch15'], d['enc'], d['rev0'] ... d['rev2'] by column
# or by row as d[0] d[1] ....

#if first file
first_ramp_i, = np.where(data_on_disk['enc'] < 20)
#skip first samples
d = data_on_disk[first_ramp_i[0]:]

# I do not want to modify the file on disk so I make a copy
# loading it into memory
outdtype=np.dtype( [('rev',np.float64),('enc',np.uint16)] + [(ch,np.float64) for ch in channels_labels] )
data = np.zeros(len(d), dtype=outdtype)

# -10 +10 V calibration
for ch in channels_labels:
    data[ch] = d[ch].astype(np.float) * 20 / 2**16 - 10 
# rev number
data['rev'] = d['rev0'].astype(np.float) + 256 * d['rev1'] + 256**2 * d['rev2']
data['enc'] = d['enc']

# removing revolutions with not all samples
start_of_revs, = np.where(data['enc']<ENC_START_TRIGGER)
samples_per_rev = np.diff(start_of_revs)
invalid_revs, = np.where(samples_per_rev != SEC_PER_REV)

for i in invalid_revs[::-1]:
    data = np.concatenate([data[:start_of_revs[i]],data[start_of_revs[i+1]:]])

# splitting data as list of revolutions
start_of_revs, = np.where(data['enc']<ENC_START_TRIGGER)
data_splitted = np.split(data[:start_of_revs[-1]], start_of_revs[1:-1])

# defining demodulated dtype
ch_dtype = np.dtype( [('T',np.float),('Q',np.float),('U',np.float)] )
demod_dtype = np.dtype( [('rev',np.float)] + [(ch,ch_dtype) for ch in channels_labels] )
demod = np.zeros(len(data_splitted), dtype=demod_dtype)

def square_wave(total_points, period, phase=0, U=False):
    '''Square wave [+1,-1]''' 
    eighth = math.floor(total_points/period)
    if U:
        phase += eighth/2
    commutator = np.array([])
    for i in range(period):
        sign = 1
        if i % 2:
            sign = -1
        commutator = np.concatenate([commutator, sign * np.ones(eighth)])
    return np.roll(commutator,int(phase))

for rev,rev_array in enumerate(data_splitted):
    for ch in channels_labels:
        channel_phase = 0
        q_commutator = square_wave(SEC_PER_REV, period=8, phase=channel_phase)
        u_commutator = square_wave(SEC_PER_REV, period=8, phase=channel_phase, U=True)

        demod['rev'][rev] = rev_array['rev'][0]
        demod[ch]['T'][rev] = np.mean(rev_array[ch])
        demod[ch]['Q'][rev] = np.mean(rev_array[ch] * q_commutator)
        demod[ch]['U'][rev] = np.mean(rev_array[ch] * u_commutator)

# plotting encoder
#plt.plot(d['enc'][:1000],label='encoder')
#plt.legend()
#plt.grid()
#plt.xlabel('Sample number')
#plt.show()
