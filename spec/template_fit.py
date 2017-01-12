# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:49:38 2016

@author: yaoyuhan
"""
import numpy as np
import collections
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.io import fits
from numpy.polynomial import polynomial as P
from scipy.optimize import leastsq
from starpy.spec.read_spectrum import read_spectrum


def _pixel_shift(wv_tg, wv_tp, fl_tp):
    """
    Return the wave and flux of template spectrum at target's pixel.
    """
    f = interpolate.interp1d(wv_tp, fl_tp)
    fl_tp_shifted = f(wv_tg)
    return fl_tp_shifted
   
   
def _coeff(flux_tp, flux_tg, deg=2):
    c = P.polyfit(flux_tp, flux_tg, deg)
    return c


def _poly_fit(wave, flux_tp, flux_tg, c):
    n = len(c)
    flux_tp_fitted = np.zeros(flux_tp.shape)
    for j in xrange(n):
        flux_tp_fitted = flux_tp_fitted + c[j]*(flux_tp**(j))
    return flux_tp_fitted


def _find_const(flux_tp, c):
    guess_const = np.mean(flux_tp)
    optimize_func = lambda x: flux_tp - x
    est_const = leastsq(optimize_func, guess_const)[0]
    best_const = 0.
    n = len(c)
    for j in xrange(n):
        best_const = best_const + c[j]*(est_const**(j))    
    return best_const


def _find_goodness(flux_tg_tofit, flux_fitted_shoulder):
    Q = np.std(flux_tg_tofit-flux_fitted_shoulder)
    return Q


def subtract_template(wv_tg, fl_tg, wv_tp, fl_tp, region_info):
    '''
    Parameters
    ----------
    wv_tg, fl_rg: wave and flux of the spectrum of target.
    wv_tp, fl_tp: wave and flux of the spectrum of template. 
        Target and template may not be in the same length(pixel).
        Flux of template must be predicted at target's pixel,
        since we may need to measure line indice after that.
    region_info: dict
        information about spectral line, eg:
        region_info = {'region_center':  6563.,
                       'region_left':    (6400., 6555.),
                       'region_right':   (6575., 6700.)}
        Left and right are to be poly_fitted, 
        whereas the whole region are to be predicted.
    Return
    ------
    wv_tg, fl_after_subtraction
    '''
    # Step 0 
    # It is necessary that the wave range of target is included in the wave 
    #   range of template, so we must ensure this precondition.
    wv_tg = np.array(wv_tg)
    fl_tg = np.array(fl_tg)
    wv_tp = np.array(wv_tp)
    fl_tp = np.array(fl_tp)
    ind_pre = np.all([wv_tg > wv_tp[0],
                      wv_tg < wv_tp[-1]],axis=0)
    wv_tg = wv_tg[ind_pre]
    fl_tg = fl_tg[ind_pre]
    
    # Step 1 : Preparation
    region_left = region_info['region_left']
    region_right = region_info['region_right']
    wave = wv_tg
    flux_tg = fl_tg
    flux_tp = _pixel_shift(wave, wv_tp, fl_tp)
    
    # Step 2
    # Only do polyfitting on the shoulder to derive cofficients
    ind_shoulder = np.any([np.all([wave > region_left[0],
                                   wave < region_left[1]], axis=0),
                           np.all([wave > region_right[0],
                                   wave < region_right[1]], axis=0)], axis=0)
    wave_tofit = wave[ind_shoulder]
    flux_tg_tofit = flux_tg[ind_shoulder]
    flux_tp_tofit = flux_tp[ind_shoulder]
    c = _coeff(flux_tp_tofit, flux_tg_tofit, deg=2) 
    
    # Step 3 : Subtract template and add a constant 
    # 3.1 find the best fitted flux
    ind_region = np.all([wave > region_left[0],
                         wave < region_right[1]],axis=0)
    # the region is where flux are to be subtracted
    wave_region = wave[ind_region]
    flux_tg_region = flux_tg[ind_region]
    flux_tp_region = flux_tp[ind_region]
    flux_fitted = _poly_fit(wave_region, flux_tp_region, 
                            flux_tg_region, c)
    flux_fitted_shoulder = _poly_fit(wave_tofit, flux_tp_tofit,
                                     flux_tg_tofit, c)
    Q = _find_goodness(flux_tg_tofit, flux_fitted_shoulder)

    # 3.2 find the constant
    const = _find_const(flux_tp_region, c)
    # and finally...
    flux_residual = flux_tg_region - flux_fitted
    flux_ulti = flux_residual + const
    return collections.OrderedDict([
            ('goodness',       Q),
            ('wave_region',    wave_region),
            ('flux_est',       flux_ulti),
            ('flux_tg_region', flux_tg_region),
            ('flux_fitted',    flux_fitted)])
            
            
def best_template_subtraction(wv_tg, fl_tg, wv_tp_list, fl_tp_list, region_info,
                              star_type = 1):
    '''
    If carbon star: set star_type = 0.
    '''
    wv_tp_list = np.array(wv_tp_list)
    fl_tp_list = np.array(fl_tp_list)
    n_tp = len(wv_tp_list)
    Q = []
    for i in range(n_tp):
        wv_tp = wv_tp_list[i]
        fl_tp = fl_tp_list[i]
        q = subtract_template(wv_tg, fl_tg, wv_tp, fl_tp, region_info)['goodness']
        Q.append(q) 
    Q = np.array(Q)
    best_ind = np.argsort(Q)[0]
    wv_tp = wv_tp_list[best_ind]
    fl_tp = fl_tp_list[best_ind]
    info = subtract_template(wv_tg, fl_tg, wv_tp, fl_tp, region_info)
    if star_type == 0:
        best_ind = 7
    info.update([('template_used', best_ind),
                 ('center',        region_info['region_center'])])
    return info


def MCgiant_tp_list_init():
    wv_tp_list_m = fits.open('/Users/yaoyuhan/Documents/data/template/wave.fits')[0].data
    fl_tp_list_m = fits.open('/Users/yaoyuhan/Documents/data/template/flux.fits')[0].data    
    n_tp_pixel = len(wv_tp_list_m[0])
    wv_tp_c = wv_tp_list_m[0]
    fl_tp_c = np.ones(n_tp_pixel)
    wv_tp_list = np.vstack([wv_tp_list_m, wv_tp_c])
    fl_tp_list = np.vstack([fl_tp_list_m, fl_tp_c])
    return wv_tp_list, fl_tp_list
    
    
if __name__ == '__main__':
    print('')
    print('@Yayaha: start to test the module ...')
    wv_tp_list, fl_tp_list = MCgiant_tp_list_init()
    region_info_Halpha = {'region_center':  6563.,
                          'region_left':    (6400., 6555.),
                          'region_right':   (6575., 6700.)}   
    fp = '/Users/yaoyuhan/Documents/data/Lamo_all_Mira/spec-57009-GAC058N16B1_sp12-181.fits'
    wv_tg = read_spectrum(fp)['wave']
    fl_tg = read_spectrum(fp)['flux'] 
    info = best_template_subtraction(wv_tg, fl_tg, wv_tp_list, fl_tp_list, 
                                     region_info=region_info_Halpha,
                                     star_type = 1)   
    print fp 
    print info.keys()                                
    plt.plot(info['wave_region'],info['flux_tg_region'])
    plt.plot(info['wave_region'],info['flux_fitted'])
    plt.plot(info['wave_region'],info['flux_est'])
    print('@Yayaha: OK')