import numpy as np
from ConfigParser import ConfigParser

config = dict(  NCHAN=16, 
                SEC_PER_REV=256, 
                ENC_START_TRIGGER=15
             )


# Names of the data channels:
channels_labels = ['ch%d' % i for i in range(config['NCHAN'])]

# Read phases
phases = ConfigParser(dict(((ch,'0') for ch in channels_labels)))
phases.read('phases.cfg')

# Structure of the data we read from the .dat files:
dat_dtype = np.dtype( [(ch,np.uint16) for ch in self.channels_labels] + 
                      [('enc',np.uint16)]+[('dummy',np.uint16)]  + 
                      [('rev%d' % i,np.uint16) for i in range(3)])

# Structure of revdata before demodulation
# in volts
rev_dtype = np.dtype( [('rev',np.long)] + 
                      [(ch,np.float,self.config['SEC_PER_REV']) for ch in self.channels_labels] )
# in ADU
rev_dtype_adu = np.dtype( [('rev',np.long)] + 
                      [(ch,np.uint16,self.config['SEC_PER_REV']) for ch in self.channels_labels] )

# Structure of a single channel
ch_dtype = np.dtype( [('T',np.float),('Q',np.float),('U',np.float)] )

# Structure of the output demodulated data
demod_dtype = np.dtype( [('rev',np.float)] + [(ch,ch_dtype) for ch in channels_labels] )

