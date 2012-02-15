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
