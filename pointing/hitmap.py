import numpy as np
import pyfits
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

if __name__ == '__main__':

    freq = 10
    f=pyfits.open('data/eq_pointing_%d.fits' % freq)
    ch = 'CHANNEL_1'
    p=f[ch].data
    for NSIDE in [64, 128, 256] :
        pix=hp.ang2pix(NSIDE,p['THETA'],p['PHI'])
        hits=hp.ma(pix2map(pix, NSIDE))
        hits.mask = hits==0
        hp.mollzoom(hits.filled(),min=0,max=int(hits.mean()*2), title="Hitmap %dGHz %s NSIDE %d" % (freq,ch, NSIDE))
