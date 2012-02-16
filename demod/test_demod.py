import unittest
import numpy as np

from .demod import *

def simple_demod(revolution_array, number_of_phases):
    revolution_length = len(revolution_array)
    phase_length = revolution_length/number_of_phases
    out = 0
    sign = +1
    for i in range(0, revolution_length, phase_length):
        out += sign * np.mean(revolution_array[i:i+phase_length])
        sign *= -1
    return out/number_of_phases

class TestDatParsing(unittest.TestCase):
    
    def setUp(self):
        """Creates a fake dataset with 2 good revolutions"""
        self.revdata = np.zeros(2, dtype=rev_dtype)
        self.revdata['ch15'][0] =  np.arange(256)
        self.revdata['ch15'][1] =  np.arange(256) + 10


    def test_demodulate(self):

        demod_data = demodulate(self.revdata, 10, number_of_phases=8)
        self.assertEqual(len(demod_data), 2)
        self.assertEqual(demod_data.dtype, demod_dtype)
        self.assertEqual(demod_data['ch15']['Q'][0], simple_demod(self.revdata['ch15'][0], 8))
        self.assertEqual(demod_data['ch15']['Q'][1], simple_demod(self.revdata['ch15'][1], 8))

