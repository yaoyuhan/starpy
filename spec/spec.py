# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:49:38 2016

@author: yaoyuhan
"""
from starpy.spec.read_spectrum import read_spectrum
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

    
class Spec(Table):

    def __init__(self, *args, **kwargs):
        super(Spec, self).__init__(*args, **kwargs)

    def norm_spec_pixel(self, norm_wave):
        sub_nearest_pixel = np.argsort(np.abs(self['wave']-norm_wave))[0]
        self['flux'] /= self['flux'][sub_nearest_pixel]
        return self

    def norm_spec_median(self):
        self['flux'] /= np.median(self['flux'])
        return self
        
    def norm_spec_max(self):
        self['flux'] /= max(abs(np.max(self['flux'])),abs(np.min(self['flux'])))
        return self

    def extract_chunk_wave_interval(self, wave_intervals=None):
        """ return spec chunk in a given wavelength interval """
        if wave_intervals is None:
            return self
        else:
            spec_chunks = []
            # should use ind (not sub)
            for wave_interval in wave_intervals:
                ind_wave_interval = np.logical_and(
                    self['wave'] >= wave_interval[0],
                    self['wave'] <= wave_interval[1])
                spec_chunks.append(self[ind_wave_interval])
            return spec_chunks
    
    def shift_to_wavenumber(self):
        self['wave'] = 1./self['wave']*10**8
        return self
            

if __name__ == '__main__':
    fp = '/Users/yaoyuhan/Documents/data/Lamo_all_Mira/spec-56958-M31017N47B1_sp02-105.fits'
    spec = Spec(read_spectrum(fp))
    wave_intervals=np.array([[4000,5000],[6000,7000]])
    spec=spec.extract_chunk_wave_interval(wave_intervals)[1]
    spec.norm_spec_pixel(6100)
    plt.plot(spec['wave'],spec['flux'])