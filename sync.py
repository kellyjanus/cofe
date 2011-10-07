#!/usr/bin/env python
import sys
import os.path
from optparse import OptionParser

from groundstation.reshape_fits import reshape_folder, concatenate_fits
from demod import demod_concatenate_fits 
from synclib import ServoSciSync

usage = """usage: %prog [options] day1, day2

Level 1 COFE sync, processes servo and scientific data to produce synchronized Level1 data
"""

parser = OptionParser(usage=usage)
parser.add_option("-f", "--folder", dest="folder",
                  help="Data base folder, expects 10_GHz_Data/ 15_GHz_Data/ and Servo_Test_Data/ subfolders", metavar="DIR")
parser.add_option("-d", "--demod", dest="demod",
                  help="Force demodulation of scientific data even if demodulated data are already available", action="store_true",default=False)
(options, args) = parser.parse_args()


if len(args) < 1:
    parser.print_help()
    sys.exit(1)

for day in args:
    print(">>>>>>>>>>>>>>>>> DAY %s" % day)

    print("PROCESSING SERVO DATA")
    folder = os.path.join(options.folder, "servo", day)
    print(folder)
    reshape_folder(folder)
    concatenate_fits(folder)

    print("PROCESSING SCIENTIFIC DATA")

    freqs = []
    for freq in [10, 15]:

        folder = os.path.join(options.folder, "%d" % freq, day)
        if os.path.exists(folder):
            freqs.append(freq)
            if (not os.path.exists(folder + '.fits')) or options.demod:
                demod_concatenate_fits(folder)

    if freqs:
        print("SYNCHRONIZING")
        s = ServoSciSync(base_folder = options.folder, day = day, freq = freqs)
        s.run()
    else:
        print("Scientific data folders not found")
