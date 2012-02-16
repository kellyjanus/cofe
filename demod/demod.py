from __future__ import division
import os.path
import numpy as np
import math
import pyfits
import logging as l
from exceptions import Exception
from optparse import OptionParser
from glob import glob

from . import utils
from . import datparsing
from .dtypes import *

"""This file provides functions to perform the demodulation of the
telescopes' .dat files."""

# Here is a lengthy explanation of what the demodulation process is and
#  why we do it that way:
# Okay. So, we have this thing called a half-wave plate, which is a metal
#  plate with a bunch of parallel wires across the front. Its effect is
#  this: light comes in. Part of it is unpolarized, and that is reflected
#  off the plate normally. Part of it is completely polarized, and the
#  direction of polarization is mirrored around the plane (a) perpendicular
#  to the disk and (b) including the wire where the light is hitting.
# Because of how mirroring works, this means that the direction of the
#  polarization we're measuring rotates twice for every time the half-wave
#  plate rotates. On top of that, since light polarized "leftwards"
#  is the same as light polarized "rightwards", we get four cycles of
#  up-down-ness to left-right-ness.
# Now, from that information, we want to calculate two numbers, Q and U.
#  Q is (vertical polarization) minus (horizontal polarization),
#  so in one revolution it goes (high,low,high,low,high,low,high,low):
#  one (high,low) for each time the measured polarization direction
#  goes from vertical to horizontal.
#  U is (up-rightish polarization) minus (down-rightish polarization),
#  so it's KIND OF like Q phase-shifted forward by half a (high,low)
#  transition.
# Now, to calculate the Q and U of the light BEFORE it gets weirded up
#  by the half-wave plate, we multiply Q by a [+1,-1] square wave to
#  approximately fix the rotation that the half-wave plate puts on,
#  the square wave having the period of the (high,low) transitions.
#  The square wave for U is the same as for Q, just phase-shifted by
#  half a hump.
# We take the mean of (data * square wave) to approximate the Q and U.
#  We also want to find T, intensity, but that's easy: just the mean
#  value of the data channel over the revolution.
#
# So! In total, what we do is this:
#  Read revolutions from a .dat file.
#  Make two copies of each channel's data.
#  Multiply each copy by a square wave.
#  Average the multiplied value over each revolution (and average the
#    unmultiplied original value).
#  Return an array 49 elements wide:
#    (revolution number + mean TQU of each channel) for each revolution.

def demodulate(data, freq, number_of_phases=8):
    """Demodulates input revdata dataset created by datparsing.create_revdata

    Parameters
    ----------
    data : ndarray
        input data of dtype rev_dtype
    freq : int
        10 or 15, required to get phases
    number_of_phases : int, optional
        number of phases to divide the revolution into
        typically 8 to get Q and U

    Returns
    -------
    demod_data : ndarray
        demodulated data of dtype demod_dtype
    """
    demod_data = np.zeros(len(data), dtype=demod_dtype)
    demod_data['rev'] = data['rev']
    for ch in channels_labels:
        calibdata = data[ch]
        channel_phase = phases.getint('%dGHz' % freq, ch)
        q_commutator = square_wave(config['SEC_PER_REV'], number_of_phases, phase=channel_phase)
        u_commutator = square_wave(config['SEC_PER_REV'], number_of_phases, phase=channel_phase, U=True)
        demod_data[ch]['T'] = np.mean(calibdata,axis=1)
        demod_data[ch]['Q'] = np.mean(calibdata*q_commutator,axis=1)
        demod_data[ch]['U'] = np.mean(calibdata*u_commutator,axis=1)
    return demod_data

def demodulate_dat(filename, freq):
    """Reads, reshapes and demodulate .dat file

    Parameters
    ----------
    filename : str
        .dat filename
    freq : int
        10 or 15, required to get phases

    Returns
    -------
    demod_data : ndarray
        demodulated data of dtype demod_dtype
    """
    return demodulate(datparsing.create_revdata(datparsing.open_raw(filename)), freq=freq)

def write_fits(demod_data, outfilename):
    """Write the demod data dictionary or compound array to a fits file
    
    Parameters
    ----------
    demod_data : compound array 
        data with demod_dtype datatype
    outfilename : str
        output filename
    """
    print('Writing %s to disk' % outfilename)
    tqu = ['T','Q','U']
    cols = []
    cols.append(pyfits.Column(name='rev', format='E', array=demod_data['rev']))
    for ch in channels_labels:
        for t in tqu:
            cols.append(pyfits.Column(name="%s_%s" % (ch,t), format='E', array=demod_data[ch][t]))
    hdu = pyfits.new_table(cols)
    hdu.header.update('EXTNAME',  'DATA')
    hdu.writeto(outfilename, clobber=True)
    print('%s written to disk' % outfilename)

if __name__ == '__main__':
    usage = '''usage: %prog [options] file_or_folder
    
    Demodulate a list of .dat files'''
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False)
    parser.add_option("-f", "--freq", action="store", type="int", dest='freq')
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        exit(1)
    if options.verbose:
        l.basicConfig(level=l.DEBUG)
    else:
        l.basicConfig(level=l.WARNING)
    for f in args:
        demod_data = demodulate_dat(f, options.freq)
        write_fits(demod_data, f.replace('.dat', '.fits'))
