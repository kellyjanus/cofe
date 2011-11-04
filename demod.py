from __future__ import division
import os.path
import numpy as np
import math
import pyfits
import logging as l
from exceptions import Exception
from glob import glob

channels_labels = ['ch%d' % i for i in range(16)]
ch_dtype = np.dtype( [('T',np.float),('Q',np.float),('U',np.float)] )
demod_dtype = np.dtype( [('rev',np.float)] + [(ch,ch_dtype) for ch in channels_labels] )
l.basicConfig(level = l.DEBUG)

class InvalidFileException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def write_fits(demod_data, outfilename):
    print('Writing %s to disk' % outfilename)
    tqu = ['T','Q','U']
    cols = []
    cols.append(pyfits.Column(name='rev', format='E', array=demod_data['rev']))
    for ch in channels_labels:
        for t in tqu:
            cols.append(pyfits.Column(name="%s_%s" % (ch,t), format='E', array=demod_data[ch][t]))
    hdu = pyfits.new_table(cols)
    hdu.header.update('EXTNAME',  'DATA')
    hdu.writeto(outfilename, clobber=True)
    print('%s written to disk' % outfilename)

def square_wave(total_points, period, phase=0, U=False):
    '''Square wave [+1,-1]''' 
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

def read_raw_data(filenames=None):
    if filenames is None:
        import tkFileDialog
        filenames = tkFileDialog.askopenfilenames()
    dr = RawDataReader(filenames)
    dr.concat_files()
    return dr.data

def concatenate_fits(filenames):

    demod_data = np.zeros(0, dtype=demod_dtype)
    outfilename = "%s-%s.fits" % (os.path.basename(filenames[0]).split('.')[0], os.path.basename(filenames[1]).split('.')[0])
    for f in filenames:
        ff = pyfits.open(f)
        file_demod = np.zeros(len(ff[1].data), dtype=demod_dtype)
        for ext in ff[1:]:
            print(ext.name)
            if (ext.name.startswith('CH')):
                for pol in ['T','Q','U']:
                    file_demod[ext.name.lower()][pol] = ext.data.field(pol)
            else:
                file_demod[ext.name.lower()] = ext.data.field(0)
        demod_data = np.concatenate([demod_data, file_demod])
    write_fits(demod_data, outfilename)

def open_demod_fits(f):
    ff = pyfits.open(f)
    file_demod = np.zeros(len(ff[1].data), dtype=demod_dtype)
    for col in ff[1].data.dtype.names:
        if (col.startswith('ch')):
            ch, pol = col.split('_')
            file_demod[ch][pol] = ff[1].data.field(col)
        else:
            file_demod[col] = ff[1].data.field(0)
    ff.close()
    return file_demod

def demod_concatenate_fits(folder):

    demod_data = np.zeros(0, dtype=demod_dtype)
    outfilename = "%s.fits" % (folder)
    
    files = sorted(glob(os.path.join(folder, '*.dat')))
    n_files = len(files)
    for n,f in enumerate(files):
        print('************** processing file %d/%d' % (n+1, n_files))
        fits_file = f.replace('.dat','.fits')
        if os.path.exists(fits_file):
            print('Fits file already converted')
            demod_data = np.concatenate([demod_data, open_demod_fits(fits_file)])
        else:
            try:
                d = DataDemod(f)
                d.run(write=True)
                demod_data = np.concatenate([demod_data, d.demod])
            except InvalidFileException as e:
                print('****************SKIPPING file %s: %s' % (f, e.value))
            
    write_fits(demod_data, outfilename)

class RawDataReader(object):

    config = dict(  NCHAN=16, 
                    SEC_PER_REV=256, 
                    ENC_START_TRIGGER=15
                 )

    def __init__(self, filenames):
        self.filenames = sorted(filenames)

    @property
    def channels_labels(self):
        return ['ch%d' % i for i in range(self.config['NCHAN'])]

    def load_raw(self, filename):
        l.info('Loading raw file %s' % filename)
        cofedtype=np.dtype( [(ch,np.uint16) for ch in self.channels_labels] + [('enc',np.uint16)]+[('dummy',np.uint16)]  + [('rev%d' % i,np.uint16) for i in range(3)])
        return np.memmap(filename,dtype=cofedtype,mode='r')

    def create_dataset(self, raw_data):
        '''Deletes invalid revolutions and shapes the array on revolutions'''
        l.info('Checking for invalid revolutions')
        start_of_revs, = np.where(raw_data['enc']<self.config['ENC_START_TRIGGER'])
        d = np.array(raw_data[start_of_revs[0]:start_of_revs[-1]].copy())
        start_of_revs, = np.where(d['enc']<self.config['ENC_START_TRIGGER'])
        samples_per_rev = np.diff(start_of_revs)
        invalid_revs, = np.where(samples_per_rev != self.config['SEC_PER_REV'])
        if len(invalid_revs) > 0:
            l.info('Removing invalid revolutions (index from beginning of file): %s' % invalid_revs)
        else:
            l.info('No invalid revolutions')
        for i in invalid_revs[::-1]:
            #d = np.concatenate([d[:start_of_revs[i]],d[start_of_revs[i+1]:]])
            d = np.delete(d, np.s_[start_of_revs[i]:start_of_revs[i+1]])
        start_of_revs, = np.where(d['enc']<self.config['ENC_START_TRIGGER'])
        samples_per_rev = np.diff(start_of_revs)
        invalid_revs, = np.where(samples_per_rev != self.config['SEC_PER_REV'])
        outdtype=np.dtype( [('rev',np.float64)] + [(ch,np.uint16,self.config['SEC_PER_REV']) for ch in self.channels_labels] )
        if len(invalid_revs) == len(start_of_revs):
            l.error('NO VALID DATA IN FILE')
            data = np.zeros(0, dtype=outdtype)
        else:
            data = np.zeros(len(d)/self.config['SEC_PER_REV'], dtype=outdtype)
            d_rev = d[::self.config['SEC_PER_REV']]
            data['rev'] =  d_rev['rev0'].astype(np.float64) + self.config['SEC_PER_REV'] * d_rev['rev1'] + self.config['SEC_PER_REV']**2 * d_rev['rev2']
            for ch in self.channels_labels:
                data[ch] = d[ch].reshape((-1, self.config['SEC_PER_REV']))
        return data

    def process_raw(self, filename):
        return self.create_dataset(self.load_raw(filename))

    def concat_files(self):
        self.data = np.concatenate(map(self.process_raw, self.filenames))

class DataDemod(object):

    config = dict(  NCHAN=16, 
                    SEC_PER_REV=256, 
                    ENC_START_TRIGGER=15
                 )

    def __init__(self, filename):
        self.filename = filename
        self.outfilename = self.filename.replace('.dat','.fits')

    def create_dataset(self):
        '''Deletes invalid revolutions and shapes the array on revolutions'''
        l.info('Checking for invalid revolutions')
        start_of_revs, = np.where(self.data_on_disk['enc']<self.config['ENC_START_TRIGGER'])
        d = np.array(self.data_on_disk[start_of_revs[0]:start_of_revs[-1]].copy())
        start_of_revs, = np.where(d['enc']<self.config['ENC_START_TRIGGER'])
        samples_per_rev = np.diff(start_of_revs)
        invalid_revs, = np.where(samples_per_rev != self.config['SEC_PER_REV'])
        if len(invalid_revs) > 0:
            l.info('Removing invalid revolutions (index from beginning of file): %s' % invalid_revs)
        else:
            l.info('No invalid revolutions')
        for i in invalid_revs[::-1]:
            #d = np.concatenate([d[:start_of_revs[i]],d[start_of_revs[i+1]:]])
            d = np.delete(d, np.s_[start_of_revs[i]:start_of_revs[i+1]])
        start_of_revs, = np.where(d['enc']<self.config['ENC_START_TRIGGER'])
        samples_per_rev = np.diff(start_of_revs)
        invalid_revs, = np.where(samples_per_rev != self.config['SEC_PER_REV'])
        if len(invalid_revs) == len(samples_per_rev):
            l.error('No valid data in file %s' % self.filename)
            raise InvalidFileException('No valid data')

        #data_splitted = np.delete(np.split(data_on_disk[:start_of_revs[-1]], start_of_revs[1:-1]), invalid_revs)

        outdtype=np.dtype( [('rev',np.long)] + [(ch,np.uint16,self.config['SEC_PER_REV']) for ch in self.channels_labels] )
        self.data = np.zeros(len(d)/self.config['SEC_PER_REV'], dtype=outdtype)
        d = d[:len(self.data)*self.config['SEC_PER_REV']]
        d_rev = d[::self.config['SEC_PER_REV']]
        self.data['rev'] =  d_rev['rev0'].astype(np.long) + self.config['SEC_PER_REV'] * d_rev['rev1'].astype(np.long) + self.config['SEC_PER_REV']**2 * d_rev['rev2'].astype(np.long)
        for ch in self.channels_labels:
            self.data[ch] = d[ch].reshape((-1, self.config['SEC_PER_REV']))

    @property
    def channels_labels(self):
        return ['ch%d' % i for i in range(self.config['NCHAN'])]

    def load_raw(self):
        l.info('Loading raw file')
        self.cofedtype=np.dtype( [(ch,np.uint16) for ch in self.channels_labels] + [('enc',np.uint16)]+[('dummy',np.uint16)]  + [('rev%d' % i,np.uint16) for i in range(3)])
        self.data_on_disk=np.memmap(self.filename,dtype=self.cofedtype,mode='r')

    def demodulate(self):
        l.info('Demodulating')
        self.demod = np.zeros(len(self.data), dtype=demod_dtype)
        self.demod['rev'] = self.data['rev']
        for ch in self.channels_labels:
            calibdata = self.data[ch].astype(np.float) * 20 / 2**16 - 10 
            channel_phase = 0
            q_commutator = square_wave(self.config['SEC_PER_REV'], period=8, phase=channel_phase)
            u_commutator = square_wave(self.config['SEC_PER_REV'], period=8, phase=channel_phase, U=True)
            self.demod[ch]['T'] = np.mean(calibdata,axis=1)
            self.demod[ch]['Q'] = np.mean(calibdata*q_commutator,axis=1)
            self.demod[ch]['U'] = np.mean(calibdata*u_commutator,axis=1)

    def write(self):
        write_fits(self.demod, self.outfilename)
        
    def run(self, write=True):
        l.info('Processing %s' % self.filename)
        if os.path.isfile(self.outfilename):
            l.warning('Demodulated file %s already created' % self.outfilename)
        else:
            self.load_raw()
            self.create_dataset()
            self.demodulate()
            if write:
                self.write()

if __name__ == '__main__':
    #l.basicConfig(level = l.DEBUG)
    #cofedata = DataDemod(sys.argv[1])
    #cofedata.run()
    print('THIS IS THE LIBRARY, RUN demod, not demod.py')
