import healpy as hp
import matplotlib.pyplot as plt
from mapmaking import gal2eq
import numpy as np

ft = hp.read_map('ftsumner.fits')
al = hp.read_map('alicesprings.fits')

ft[ft<0]=0
al[al<0]=0

NSIDE = 32
hit = hp.ud_grade(ft + al, NSIDE, power=-2)

noise = {8:308, 10:307, 15:163} #microK rt Sec
fsamp = 30. #Hz

#for freq, no in zip([8, 10, 15], noise):
freq = 8 
key = 'map'
if freq == 8:
    key = 'brian'
signal = gal2eq(hp.remove_monopole(hp.ud_grade(
hp.read_map('/smbmount/data1/COFE/Simulations/COFE_simulationPSMv163/cofe1/frequency_maps/coadded%s_%dGhz.fits' % (key,freq)), NSIDE),gal_cut=50))
if freq == 8:
    signal *= 1.e3

noise_map = noise[freq] / np.sqrt(hit * 1/fsamp) / 1e3

signoise = hp.ma(np.absolute(signal)/noise_map)
signoise.mask = hit == 0
hp.mollview(signoise, min=1,max=10.,coord='CG',title='Signal to noise at %d GHz' % freq,xsize=2000, norm='log')
plt.savefig('sigtonoise_%d.png' % freq)
