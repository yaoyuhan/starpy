# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:49:38 2016

@author: yaoyuhan
"""

import math
import collections
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import leastsq
from starpy.tools.walk_dir import walk_dir
from starpy.tools.coo import degreeTohms_dms


def _find_period(time, data, guess_period):
    guess_phase = 0
    guess_mean = np.mean(data)
    guess_std = 3*np.std(data)/(2**0.5)
    optimize_func = lambda x: x[0]*np.sin(time*2*math.pi/x[3]+x[1])+x[2]-data
    est_std, est_phase, est_mean, est_period = \
        leastsq(optimize_func, [guess_std, guess_phase, guess_mean, guess_period])[0]
    return est_std, est_phase, est_mean, est_period



def Fourier_curve_fitting(time_o, mag_o, p_guess=258.1, 
                          predict_start=55430, predict_end=57800):
    """
    First order only.
    Mira variables pulsate in fundamental mode.
    
    Parameters
    ----------
    time_o: array
        Observation time in original order, mjd_form.
    mag_o: array
        Observation magnitude in original order.
    """
    time_o = np.array(time_o)
    mag_o = np.array(mag_o)
    ind = np.argsort(time_o)
    time = time_o[ind]
    data = mag_o[ind]

    guess_period = p_guess
    
    # Define the function to optimize, we want to minimize the difference
    # between the actual data and our "guessed" parameters
    est_std, est_phase, est_mean, est_period = _find_period(time, data, guess_period)

    # Recreate the fitted curve using the optimized parameters
    t = np.linspace(predict_start,predict_end,(predict_end-predict_start)/10+1)
    data_fit = est_std*np.sin(t*2*math.pi/est_period+est_phase) + est_mean
    
    # Calculate period_err
    data_fit_at_obs = est_std*np.sin(time*2*math.pi/est_period+est_phase) + est_mean
    noise_std = np.std(data - data_fit_at_obs)
    num_refit_=100
    a = (data_fit_at_obs.reshape(1, -1).repeat(num_refit_, axis=0) +
         np.random.randn(num_refit_, data_fit_at_obs.size) * noise_std) 
    period_refit = []
    for i in range(num_refit_):
        period_refit.append(_find_period(time, a[i], guess_period)[3])
    period_refit = np.array(period_refit)
    period_err = np.std(period_refit)
    
    theta = pdm_theta(time, data, est_period)    
    
    return collections.OrderedDict([
            ('period',        est_period),
            ('period_err',    period_err),
            ('amplitude',     abs(est_std)*2),
            ('mean',          est_mean),
            ('mjd_obs',       time),
            ('mag_obs',       data),
            ('mjd_fit',       t),
            ('mag_fit',       data_fit),
            ('mag_fit_at_obs',data_fit_at_obs),
            ('naughty_degree',theta)])


def _cal_phase(lc_info, lamost_mjd):
    t = lc_info['mjd_fit']
    l = lc_info['mag_fit']
    period = lc_info['period']
    period_err = lc_info['period_err']
    mjd_range = np.array([lamost_mjd-period, lamost_mjd])
    ind_range = np.all([t > mjd_range[0],
                        t < mjd_range[1]],axis=0)
    time = t[ind_range]
    light = l[ind_range]
    ind_most_luminous = np.where(light==min(light))[0][0]
    t_most_luminous = time[ind_most_luminous]
    time_step = lamost_mjd - t_most_luminous
    phase = (time_step) / period  
    phase_err = (time_step/period**2)*period_err
    return phase, phase_err


def pdm_theta(t,x,p,N_B=10):
    """
    Reference
    ---------
    Stellingwerf B. F., 1978, ApJ, 224, 953
    """
    # Arange t & x in ascending sequence.
    x = np.array(x)
    t = np.array(t)
    ind = np.argsort(t)
    t = t[ind]
    x = x[ind] 
    N = len(x)
    V = np.var(x)*N/(N-1)
    t_0 = t[0]
    pre_phases = (t-t_0)/p
    phases = np.array([pre_phases[i]-int(pre_phases[i]) for i in xrange(N)])
    mag_groups = [[]for i in xrange(N_B)]
    dividor = 1./N_B    
    for i in xrange(N):
        mag_groups[int(phases[i]/dividor)].append(x[i])
    mag_groups=np.array(mag_groups)
    v = np.array([np.var(mag_groups[i]) for i in xrange(N_B)])
    v_bar = np.nanmean(v)    
    theta = v_bar/V
    return theta


def light_curve_kiso(fp, mjd_assign=False):
    'Deal with Individual filepath'
    
    lc_pool = np.array(walk_dir('/Users/yaoyuhan/Documents/data/lc_dat/'))
    KISOGP_info = fits.open('/Users/yaoyuhan/Documents/data/catalogue/KISOGP_Miras.fits')[1].data
    KISOGP_coo = KISOGP_info['Kiso_coo']  
    KISOGP_period = KISOGP_info['Period']
    hl = fits.open(fp)[0].header
    RA = hl['RA']
    DEC = hl['dec']
    if mjd_assign==False:
        lamost_mjd = hl['MJD']
    else:
        lamost_mjd = 57783
    
    # Difficult matching
    n_pool = len(lc_pool)
    lc_pool_ra_pre = np.array([x[38:42] for x in lc_pool])
    lc_pool_dec_pre = np.array([x[47:51] for x in lc_pool])
    ra_pre = degreeTohms_dms(RA, DEC)[0][:4]
    dec_pre = degreeTohms_dms(RA, DEC)[1][:4]
    lc_pool_coo_pre = []
    for i in xrange(n_pool):
        lc_pool_coo_pre.append(lc_pool_ra_pre[i]+lc_pool_dec_pre[i])
    lc_pool_coo_pre = np.array(lc_pool_coo_pre)
    coo_pre = ra_pre + dec_pre    
    if coo_pre not in lc_pool_coo_pre:
        return
    else:
        # ind_path is the position of the lightcurve
        # ind_coo is the position in 783 kiso miras
        ind_path = np.where(lc_pool_coo_pre == coo_pre)[0][0]
        lc_path = lc_pool[ind_path]
        kiso_coo = lc_path[38:53]
        data = Table.read(lc_path,  format = 'ascii')  
        time_obs = data['col1']
        mag_obs = data['col2']
        ind_coo = np.where(KISOGP_coo == kiso_coo)[0][0]
        period_guess = KISOGP_period[ind_coo]
        lc_info = Fourier_curve_fitting(time_obs, mag_obs, p_guess=period_guess)
        f = interpolate.interp1d(lc_info['mjd_fit'], lc_info['mag_fit'])   
        lamost_mag = float(f(lamost_mjd))
        phase, phase_err = _cal_phase(lc_info, lamost_mjd)
        lc_info.update([('lamost_mjd',          lamost_mjd),
                        ('lamost_mag',          lamost_mag),
                        ('lamost_phase',        phase),
                        ('lamost_phase_err',    phase_err),
                        ('kiso_coo',            kiso_coo),
                        ('RA',                  RA),
                        ('DEC',                 DEC)])
        return lc_info
        
        
if __name__ == '__main__':
    print('')
    print('@Yayaha: start to test the module ...')
    fp = '/Users/yaoyuhan/Documents/data/Lamo_kiso_Mira/spec-57324-NGC7788_4_sp05-143.fits'
    lc = light_curve_kiso(fp,mjd_assign=True)  
    plt.plot(lc['mjd_obs'],lc['mag_obs'],'o',color='orange')
    plt.plot(lc['mjd_fit'],lc['mag_fit'],color='b')
    print lc.keys()
    print('@Yayaha: OK')