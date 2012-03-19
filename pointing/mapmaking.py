import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pycfitsio as fits
import itertools
import pyfits
import healpy as hp

def pix2map(pix, nside, tod=None, dividebyhits=True):
    """Pixel array to hitmap, if TOD with same lenght of PIX is provided, 
    it is binned to a map"""
    #TODO test case
    pix = pix.astype(np.int)
    ids = np.bincount(pix, weights=None)
    hitmap = np.ones(hp.nside2npix(nside)+1) * hp.UNSEEN
    hitmap[:len(ids)] = ids
    hitmap = hp.ma(hitmap[:-1])
    if tod is None:
        return hitmap
    else:
        ids_binned = np.bincount(pix, weights=tod)
        binned = np.ones(hp.nside2npix(nside)+1) * hp.UNSEEN
        binned[:len(ids_binned)] = ids_binned
        binned = hp.ma(binned[:-1])
        if dividebyhits:
            binned /= hitmap
        return hitmap, binned

def weightmap(summap, invM):
    binmap = summap.copy()
    for  x in range(len(summap)):
        if invM[x][0,0]>0:
            binmap[x] = np.dot(summap[x], invM[x])
    return binmap

def get_qu_weights(pa):
    #cos2pa = np.cos(2 * pa)
    #sin2pa = np.sin(2 * pa)
    cos2pa = np.cos(pa)**2 - np.sin(pa)**2
    sin2pa = 2*np.sin(pa)*np.cos(pa)
    q_channel_w = {'Q': -1 * cos2pa, 'U': -1 * sin2pa}
    u_channel_w = {'Q': -1 * sin2pa, 'U': cos2pa}
    return q_channel_w, u_channel_w

def qubin(pix, NSIDE, q_channel, u_channel, pa):

    npix = hp.nside2npix(NSIDE)
    q_channel_w, u_channel_w = get_qu_weights(pa)

    #create M
    M = {}
    for c1,c2 in itertools.combinations_with_replacement('QU',2):
        print (c1,c2)
        M[(c1,c2)] = pix2map(pix, NSIDE, q_channel_w[c1] * q_channel_w[c2] + u_channel_w[c1] * u_channel_w[c2], dividebyhits=False)[1].filled()

    #invert M
    invM = np.array([ 
            [[ M[('Q','Q')][x], M[('Q','U')][x] ],
                           [ M[('Q','U')][x], M[('U','U')][x] ]] for x in range(npix) ])
    for x,blockM in enumerate(invM):
        if blockM[0,0] != 0:
            invM[x] = np.linalg.inv(blockM)

    summap = hp.ma(np.hstack([
        pix2map(pix, NSIDE, q_channel * q_channel_w[c] + u_channel * u_channel_w[c], dividebyhits=False)[1].filled()[:,None] for c in 'QU'
        ]))
    summap.mask = summap == 0

    return M, invM, summap, weightmap(summap, invM)

def unwrap(phase):
    i,=np.where(np.diff(phase)<0)
    i+=1
    offsets = 360. * np.repeat(np.arange(len(i)+1), np.concatenate([i, [len(phase)]]) -np.          concatenate([[0],i]))
    return phase + offsets


def phasebin(nbins, az, signal):
    """Phasebin signal for each ring in a pseudomap of nbins bins
    """
    ring_edges, = np.where(np.diff(az) < -np.pi)
    nrings = len(ring_edges)
    phasebin_edges = np.linspace(-np.pi,np.pi, nbins+1)

    pseudomap = np.zeros([nbins,nrings-1],dtype=np.float32 )
    for ring in range(nrings-1):
        az_ring=az[ring_edges[ring]:ring_edges[ring+1]]
        signal_ring=signal[ring_edges[ring]:ring_edges[ring+1]]
        pseudomap[:,ring], edges = np.histogram(az_ring, bins=phasebin_edges, weights=signal_ring)
        hits, edges = np.histogram(az_ring, bins=phasebin_edges)
        pseudomap[hits>0,ring] /= hits[hits>0]
    return pseudomap

def gal2eq(m):
    Rgal2eq = hp.Rotator(coord='CG')
    npix = len(m)
    nside = hp.npix2nside(npix)
    newpix = hp.vec2pix(nside, *Rgal2eq(hp.pix2vec(nside, np.arange(npix))))
    return m[newpix]

if __name__ == '__main__':

    freq = 10
    ch_num = 5
    mag = True
    filename = 'data/eq_pointing_%d.fits' % freq
    if mag:
        filename = filename.replace('.fits','_mag.fits')
    f=pyfits.open(filename)
    ch = 'CHANNEL_%d' % ch_num
    p=f[ch].data
    NSIDE = 128
    NPIX = hp.nside2npix(NSIDE)
    pix=hp.ang2pix(NSIDE,p['THETA'],p['PHI'])
    flagsfile = 'data/flags%dghz.fits' % freq
    datafile = 'data/all_%dGHz_data_cal.fits' % freq
    
    data_ut = fits.read(datafile, 'TIME')['UT']
    data_range = slice(data_ut.searchsorted(f['TIME'].data['UT'][0]), data_ut.searchsorted(f['TIME'].data['UT'][-1])+1)
    ch_num -= 1
    flag = fits.read(flagsfile, 'CH%d' % ch_num)['FLAGS'][data_range]
    data = fits.read(datafile, 'CH%d_' % ch_num)
    if True: #FLAG
        pix[flag > 0] = NPIX
    hits=hp.ma(pix2map(pix, NSIDE))
    hits.mask = hits==0
    #hp.mollzoom(hits.filled(),min=0,max=int(hits.mean()*2), title="Hitmap %dGHz %s NSIDE %d" % (freq,ch, NSIDE))
    print "open wmap map"
    m = fits.read_map('wmap/wmap_band_iqumap_r9_7yr_K_v4.fits', field=slice(0,2+1), nest=None)
    qu_wmap = [gal2eq(hp.ud_grade(amap, NSIDE, order_in='NESTED', order_out='RING')) for amap in m]
    qu_wmap = [np.append(amap, [0]) for amap in qu_wmap]

    q_channel_w, u_channel_w = get_qu_weights(p['PSI'])
    t_channel = qu_wmap[0][pix]
    q_channel = q_channel_w['Q'] * qu_wmap[1][pix] +  q_channel_w['U'] * qu_wmap[2][pix]
    u_channel = u_channel_w['Q'] * qu_wmap[1][pix] +  u_channel_w['U'] * qu_wmap[2][pix]
    wmapM, wiaminvM, wmap_qu_summap, wmap_qu_binned = qubin(pix, NSIDE, q_channel, u_channel, p['PSI'])
    hp.mollview(wmap_qu_binned[:,0].filled(), min=-.5, max=.5, unit='mK', title='WMAP 23GHz 7y Q', xsize=2000)
    plt.savefig('wmap_q.png')
    hp.mollview(wmap_qu_binned[:,1].filled(), min=-.5, max=.5, unit='mK', title='WMAP 23GHz 7y U', xsize=2000)
    plt.savefig('wmap_u.png')
    wmap_t_binned = pix2map(pix, NSIDE, t_channel)[1]
    hp.mollview(wmap_t_binned.filled(), min=-100, max=100, unit='mK', title='WMAP 23GHz 7y I', xsize=2000)
    plt.savefig('wmap_t.png')
    cofe_t_binned = pix2map(pix, NSIDE, data['T'][data_range])[1]
    M, invM,cofe_qu_summap,  cofe_qu_binned = qubin(pix, NSIDE, data['Q'][data_range], data['U'][data_range], p['PSI'])
    hp.mollview(cofe_qu_binned[:,0].filled(), min=-.5, max=.5, unit='?', title='COFE %dGHz ch %d Q' % (freq, ch_num), xsize=2000)
    plt.savefig('cofe_q_%d_ch%d.png' % (freq, ch_num))
    hp.mollview(cofe_qu_binned[:,1].filled(), min=-.5, max=.5, unit='?', title='COFE %dGHz ch %d U' % (freq, ch_num), xsize=2000)
    plt.savefig('cofe_u_%d_ch%d.png' % (freq, ch_num))
    hp.mollview(cofe_t_binned.filled(), min=-.5, max=.5, unit='?', title='COFE %dGHz ch %d I' % (freq, ch_num), xsize=2000)
    plt.savefig('cofe_t_%d_ch%d.png' % (freq, ch_num))
    for c in 'tqu':
        subprocess.call(['convert -delay 100 '+(' wmap_%s.png ' % c)+' cofe_%s_%d_ch%d.png ' % (c, freq, ch_num)+' cofewmap_%s_%d_ch%d.gif ' % (c, freq, ch_num)],shell=True)
