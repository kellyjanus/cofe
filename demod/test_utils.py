import unittest
import numpy as np

from .utils import *

class TestUtils(unittest.TestCase):
    
    def test_squarewave(self):
        total_points = 24
        number_of_phases = 8
        out = square_wave(total_points, number_of_phases, phase=0, U=False)
        ones = np.ones(3)
        np.testing.assert_array_equal(out, np.concatenate(
            [ones, -1*ones, ones, -1*ones, ones, -1*ones, ones, -1*ones]))

    def test_squarewave_U(self):
        total_points = 32
        number_of_phases = 8
        out = square_wave(total_points, number_of_phases, phase=0, U=True)
        ones = np.ones(4)
        np.testing.assert_array_equal(out, np.concatenate(
            [[-1,-1], ones, -1*ones, ones, -1*ones, ones, -1*ones, ones, [-1,-1]]))

    def test_squarewave_phase(self):
        total_points = 32
        number_of_phases = 8
        out = square_wave(total_points, number_of_phases, phase=1, U=False)
        ones = np.ones(4)
        np.testing.assert_array_equal(out, np.concatenate(
            [[-1], ones, -1*ones, ones, -1*ones, ones, -1*ones, ones, [-1,-1,-1]]))
