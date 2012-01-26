import numpy as np
import healpy as hp

def pix2map(pix, nside, tod=None):
    """Pixel array to hitmap, if TOD with same lenght of PIX is provided, 
    it is binned to a map"""
    #TODO test case
    pix = pix.astype(np.int)
    ids = np.bincount(pix, weights=None)
    hitmap = np.ones(hp.nside2npix(nside)) * hp.UNSEEN
    hitmap[:len(ids)] = ids
    hitmap = hp.ma(hitmap)
    if tod is None:
        return hitmap
    else:
        ids_binned = np.bincount(pix, weights=tod)
        binned = np.ones(hp.nside2npix(nside)) * hp.UNSEEN
        binned[:len(ids_binned)] = ids_binned
        binned = hp.ma(binned)/hitmap
        return hitmap, binned
