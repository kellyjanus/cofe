import numpy as np
import math

def adu2volts(adu):
    """Convert ADU output into volts

    Parameters
    ----------
    adu : integer
        input ADU

    Returns
    -------
    volts : float
        value converted into volts

    Examples
    --------
    >>> adu2volts(2.**12)
    -8.75
    """
    return adu * 20./2**16 - 10

def square_wave(total_points, number_of_phases, phase=0, U=False):
    """Square wave [+1,-1]
    
    Parameters
    ----------
    total_points : int
        length of the output array
    number_of_phases : int
        number of phases to divide the total points into
    phase : int, optional
        phase in number of samples
    U : bool, optional
        if U is true, phase is shifted of half the period
        
    Returns
    -------
    out : ndarray
        output array with square wave
    """ 

    period_length = total_points // number_of_phases

    # U is shifted by half period
    if U:
        phase += period_length // 2

    # create +1, -1, +1, -1 array with 1 element for each state
    states_multiplier = np.ones(number_of_phases)
    states_multiplier[1::2] *= -1
    
    # repeat by the length of the period
    commutator = np.repeat(states_multiplier, period_length)

    # apply the phases
    out = np.roll(commutator,int(phase))
    return out
