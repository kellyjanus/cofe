import numpy as np

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
    """
    return adu * 20./2**16 - 10

def square_wave(total_points, period, phase=0, U=False):
    """Square wave [+1,-1]
    
    Parameters
    ----------
    total_points : int
        length of the output array
    period : int
        period length in number of samples
    phase : int, optional
        phase in number of samples
    U : bool, optional
        if U is true, phase is shifted of half the period
        
    Returns
    -------
    out : ndarray
        output array with square wave
    """ 

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

