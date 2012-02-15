"""This file contains tools for parsing the data contained in the .dat files
 produced by the telescopes' data acquisition code.
The most important things in this file are:
 dat_dtype: the numpy dtype of the .dat files.
 rev_dtype: the numpy dtype of the output data
 read_raw(filename): a numpy.memmap of dtype=dat_dtype into the filename
We recognize when a new revolution starts by seeing that the encoder
 value is less than ENCODER_START_TRIGGER.
A "valid" revolution is one with exactly SAMPLES_PER_REVOLUTION samples.
 Any other is damaged in some way, and recommended to be ignored."""

import numpy as np
import logging as l

from . import utils

config = dict(  NCHAN=16, 
                SEC_PER_REV=256, 
                ENC_START_TRIGGER=15
             )

# Names of the data channels:
channels_labels = ['ch%d' % i for i in range(config['NCHAN'])]

# Structure of the data we read from the .dat files:
dat_dtype = np.dtype( [(ch,np.uint16) for ch in self.channels_labels] + 
                      [('enc',np.uint16)]+[('dummy',np.uint16)]  + 
                      [('rev%d' % i,np.uint16) for i in range(3)])

# Structure of a single channel
ch_dtype = np.dtype( [('T',np.float),('Q',np.float),('U',np.float)] )
# Structure of the output demodulated data

# in volts
rev_dtype = np.dtype( [('rev',np.long)] + 
                      [(ch,np.float,self.config['SEC_PER_REV']) for ch in self.channels_labels] )
# in ADU
rev_dtype_adu = np.dtype( [('rev',np.long)] + 
                      [(ch,np.uint16,self.config['SEC_PER_REV']) for ch in self.channels_labels] )

def open_raw(filename):
    """Reads a .dat file into a memmap

    Parameters
    ----------
    filename : str
        input .dat filename

    Returns
    -------
    out : memmap
        numpy memory map array
    """
    l.info('Loading raw file %s' % filename)
    raw_data = np.memmap( filename, dtype=dat_dtype, mode='r')
    return raw_data

def create_revdata(raw_data, volts=True):
    """Deletes invalid revolutions and shapes the array on revolutions
    
    Parameters
    ----------
    raw_data : ndarray
        input array with dtype dat_dtype

    Returns
    -------
    revdata : ndarray
        reshaped output dataset
    """

    # remove partial revolutions at the beginning and end of dataset
    start_of_revs, = np.where(raw_data['enc'] < config['ENC_START_TRIGGER'])
    d = np.array(raw_data[start_of_revs[0]:start_of_revs[-1]].copy())

    # remove revolutions with bad number of samples 
    start_of_revs, = np.where(d['enc'] < config['ENC_START_TRIGGER'])
    samples_per_rev = np.diff(start_of_revs)
    invalid_revs, = np.where(samples_per_rev != config['SEC_PER_REV'])

    if len(invalid_revs) > 0:
        l.warning('Removing invalid revolutions (index from beginning of file): %s' % invalid_revs)
    else:
        l.warning('No invalid revolutions')

    # remove the samples of the bad revolutions from the array
    for i in invalid_revs[::-1]:
        d = np.delete(d, np.s_[start_of_revs[i]:start_of_revs[i+1]])

    out_dtype = rev_dtype if volts else rev_dtype_adu

    # create the rev output array
    if len(d) == 0:
        l.error('NO VALID DATA IN FILE')
        data = np.zeros(0, dtype=out_dtype)
    else:
        data = np.zeros(len(d)/config['SEC_PER_REV'], dtype=out_dtype)
        d_rev = d[::config['SEC_PER_REV']]
        data['rev'] = d_rev['rev0'].astype(np.long) + 
                      d_rev['rev1'].astype(np.long) * self.config['SEC_PER_REV'] *  + 
                      d_rev['rev2'].astype(np.long) * self.config['SEC_PER_REV']**2
        for ch in self.channels_labels:
            chdata = d[ch].reshape((-1, self.config['SEC_PER_REV']))
            if volts:
                data[ch] = utils.adu2volts(chdata)
            else:
                data[ch] = chdata 

    return data

def read_raw(filenames, volts=True):
    """Reads a list of filenames, creates revdata dataset and concatenates them
    
    Parameters
    ----------
    filenames : list
        list of .dat filenames to read
    volts : bool, optional
        whether to convert to volts or keep ADU

    Returns
    -------
    revdata : array
        reshaped concatenated array
    """
    return np.concatenate(
                [create_revdata(open_raw(f), volts) for f in filenames]
                         )
