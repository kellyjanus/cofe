"""This file contains tools for parsing the data contained in the .dat files
produced by the telescopes' data acquisition code.
"""
import numpy as np
import logging as l

from . import utils
from .dtypes import *

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
        data['rev'] = d_rev['rev0'].astype(np.long) + \
                      d_rev['rev1'].astype(np.long) * config['SEC_PER_REV'] + \
                      d_rev['rev2'].astype(np.long) * config['SEC_PER_REV']**2
        for ch in channels_labels:
            chdata = d[ch].reshape((-1, config['SEC_PER_REV']))
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
