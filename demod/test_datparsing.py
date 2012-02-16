import unittest
import numpy as np

from .datparsing import *

class TestDatParsing(unittest.TestCase):
    
    def setUp(self):
        """Creates a fake dataset with 2 good revolutions, 1 bad revolution,
        and some random samples at the beginning and the end"""
        self.raw_data = np.zeros(256 * 3 + 20, dtype=dat_dtype)
        start = 5
        self.raw_data['enc'][:start] = 100 # bad start to remove
        end = start + 256
        self.raw_data['enc'][start:end] = np.arange(2, 4098, 16) # complete rev
        self.raw_data['ch1'][start:end] = np.arange(256)
        self.raw_data['rev0'][start] = 2
        self.raw_data['rev1'][start] = 3
        self.raw_data['rev2'][start] = 4
        self.revcheck = 2 + 3 * 256 + 4 * 256**2
        start += 256; end+=251
        self.raw_data['enc'][start:end] = np.arange(2, 4098-16*5, 16) # partial rev
        self.raw_data['enc'][end:end+256] = np.arange(10, 4106, 16) # complete rev
        self.raw_data['ch1'][end:end+256] = np.arange(256) + 2
        self.raw_data['enc'][end+256:] = np.arange(10, 10+ (len(self.raw_data) - (end+256))*16, 16)

    def test_create_revdata(self):
        out = create_revdata(self.raw_data, volts=False)

        # check number of revolutions found
        self.assertEqual(len(out), 2)
        # check correct datatype
        self.assertEqual(out.dtype, rev_dtype_adu)
        # check reconstruction of revolution counter
        self.assertEqual(out['rev'][0], self.revcheck)
        # check reshaping
        np.testing.assert_array_equal(out['ch1'][0], np.arange(256))
        np.testing.assert_array_equal(out['ch1'][1], np.arange(256)+2)

    def test_create_revdata(self):
        out = create_revdata(self.raw_data, volts=True)

        # check reshaping
        np.testing.assert_array_almost_equal(out['ch1'][0], np.arange(256)* 20./2**16 - 10
        )
        np.testing.assert_array_almost_equal(out['ch1'][1], (np.arange(256)+2) * 20./2**16 - 10
        )
