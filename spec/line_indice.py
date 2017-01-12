# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:49:38 2016

@author: yaoyuhan
"""
# sys.path must include 'Users/yaoyuhan/Documents/
import numpy as np
import collections
from lmfit.models import LinearModel, GaussianModel
from starpy.spec.read_spectrum import read_spectrum


def measure_line_index(wave,
                       flux,
                       flux_err=None,
                       mask=None,
                       z=None,
                       line_info=None,
                       Gauss = False,
                       num_refit=(100,None),
                       return_type='dict'):
    '''
    Return
    ------
    line_indx = collections.OrderedDict([
            ('SN_local_flux_std',        1. / noise_std),
            ('EW_int',                   EW_int),
            ('EW_int_err',               EW_int_err),
            ('mod_gauss_center',         ~),
            ('mod_gauss_center_err',     ~),
            ('mod_gauss_center_std',     ~
            )])         
    '''
    try:
        # 0. do some input check
        # 0.1> check line_info
        line_info_keys = line_info.keys()
        assert 'line_range' in line_info_keys
        assert 'line_shoulder_left' in line_info_keys
        assert 'line_shoulder_right' in line_info_keys
        # 0.2> check line range/shoulder in spectral range
        assert np.min(wave) <= line_info['line_shoulder_left'][0]
        assert np.max(wave) >= line_info['line_shoulder_right'][0]

        # 1. get line information
        # line_center = line_info['line_center']  # not used
        line_range = line_info['line_range']
        line_shoulder_left = line_info['line_shoulder_left']
        line_shoulder_right = line_info['line_shoulder_right']
        
        # 2. shift spectra to rest-frame
        wave = np.array(wave)
        flux = np.array(flux)
        if z is not None:
            wave /= 1. + z
        
        # 3. estimate the local continuum
        # 3.1> check if minus flux, this might not be necessary.
        ind_region = np.any([
                np.all([wave > line_shoulder_left[0],
                        wave < line_shoulder_right[1]], axis=0)], axis=0)
        flux_region = flux[ind_region]
        ind_minus = np.all([flux_region<0],axis=0)
        if True in ind_minus:
            return measure_line_index_null_result()
        
        # 3.2> shoulder wavelength range
        ind_shoulder = np.any([
            np.all([wave > line_shoulder_left[0],
                    wave < line_shoulder_left[1]], axis=0),
            np.all([wave > line_shoulder_right[0],
                    wave < line_shoulder_right[1]], axis=0)], axis=0)
        wave_shoulder = wave[ind_shoulder]
        flux_shoulder = flux[ind_shoulder]
        
        # 3.3> integrated/fitted wavelength range
        ind_range = np.logical_and(wave > line_range[0], wave < line_range[1])
        wave_range = wave[ind_range]
        flux_range = flux[ind_range]
        
        # 4. linear model
        mod_linear = LinearModel(prefix='mod_linear_')
        par_linear = mod_linear.guess(flux_shoulder, x=wave_shoulder)
        # ############################################# #
        # to see the parameter names:                   #
        # model_linear.param_names                      #
        # {'linear_fun_intercept', 'linear_fun_slope'}  #
        # ############################################# #
        out_linear = mod_linear.fit(flux_shoulder,
                                    par_linear,
                                    x=wave_shoulder,
                                    method='leastsq')
        
        # 5. estimate continuum
        cont_shoulder = out_linear.best_fit
        noise_std = np.std(flux_shoulder / cont_shoulder)
        cont_range = mod_linear.eval(out_linear.params, x=wave_range)
        resi_range = 1 - flux_range / cont_range
        
        # 6.1 Integrated EW 
        # estimate EW_int
        wave_diff = np.diff(wave_range)
        wave_step = np.mean(np.vstack([np.hstack([wave_diff[0], wave_diff]),
                                       np.hstack([wave_diff, wave_diff[-1]])]),
                            axis=0)
        EW_int = np.dot(resi_range, wave_step)
        
        # estimate EW_int_err
        num_refit_=num_refit[0]
        if num_refit_ is not None and num_refit_>0:
            EW_int_err = np.std(np.dot(
                (resi_range.reshape(1, -1).repeat(num_refit_, axis=0) +
                 np.random.randn(num_refit_, resi_range.size) * noise_std),
                wave_step))
        else:
            EW_int_err = np.nan
        
        line_indx = collections.OrderedDict([
            ('SN_local_flux_std',        1. / noise_std),
            ('EW_int',                   EW_int),
            ('EW_int_err',               EW_int_err),
            ('mod_gauss_center',         np.nan),
            ('mod_gauss_center_err',     np.nan),
            ('mod_gauss_center_std',     np.nan)])           
        
        # 7.2 Gaussian model
        if Gauss == True:
            mod_gauss = GaussianModel(prefix='mod_gauss_')
            par_gauss = mod_gauss.guess(resi_range, x=wave_range)
            out_gauss = mod_gauss.fit(resi_range, par_gauss, x=wave_range)
            mod_gauss_center = out_gauss.params[mod_gauss.prefix + 'center'].value                
            mod_gauss_center_err = out_gauss.params[mod_gauss.prefix + 'center'].stderr
            
            line_indx.update([('mod_gauss_center',        mod_gauss_center),
                              ('mod_gauss_center_err',    mod_gauss_center_err)])    
        
            num_refit_=num_refit[1]
            if num_refit_ is not None and num_refit_ > 2:
                out_gauss_refit_center = np.zeros(num_refit_)
                for i in range(num_refit_):
                    resi_range_with_noise = resi_range + \
                                    np.random.randn(resi_range.size) * noise_std
                    out_gauss_refit = mod_gauss.fit(resi_range_with_noise,
                                            par_gauss,
                                            x=wave_range)
                    out_gauss_refit_center[i]=\
                        out_gauss_refit.params[mod_gauss.prefix + 'center'].value
                mod_gauss_center_std = np.nanstd(out_gauss_refit_center)   
                line_indx.update([('mod_gauss_center_std',     mod_gauss_center_std)])
        # if necessary, convert to array
        # NOTE: for a non-ordered dict the order of keys and values may change!
        if return_type == 'array':
            return np.array(line_indx.values())

        return line_indx
    except Exception:
        return measure_line_index_null_result(return_type)


def measure_line_index_null_result(return_type):
    line_indx = collections.OrderedDict([
                  ('SN_local_flux_std',        np.nan),
                  ('EW_int',                   np.nan),
                  ('EW_int_err',               np.nan),
                  ('mod_gauss_center',         np.nan),
                  ('mod_gauss_center_err',     np.nan),
                  ('mod_gauss_center_std',     np.nan)])
    if return_type == 'array':
        return np.array(line_indx.values())
    return line_indx
    
    
def measure_Balmer_series(filepath, num_refit=(None,None), Gauss=False, return_type='array'):
    line_info_Hdelta = {'line_center':         4102,
                        'line_range':          (4096, 4109),
                        'line_shoulder_left':  (4080, 4094),
                        'line_shoulder_right': (4112, 4127)}
                        
    line_info_Hgamma = {'line_center':         4340,
                        'line_range':          (4337, 4347),
                        'line_shoulder_left':  (4320, 4335),
                        'line_shoulder_right': (4348.57, 4351.88)}
                        
    line_info_Hbeta = {'line_center':         4861,
                       'line_range':          (4857, 4868),
                       'line_shoulder_left':  (4845, 4855),
                       'line_shoulder_right': (4875, 4880)}                        
                        
    line_info_Halpha = {'line_center':         6563,
                        'line_range':          (6548.00, 6578.00),
                        'line_shoulder_left':  (6505.00, 6540.00),
                        'line_shoulder_right': (6580.00, 6610.00)}                    
           
    spec = read_spectrum(filepath)  
    wave = spec['wave']
    flux = spec['flux']   
    line_indx_Hdelta = measure_line_index(wave,flux,line_info=line_info_Hdelta,
                                   num_refit = num_refit, Gauss=Gauss,
                                   return_type=return_type)
                                   
    line_indx_Hgamma = measure_line_index(wave,flux,line_info=line_info_Hgamma,
                                   num_refit = num_refit, Gauss=Gauss,
                                   return_type=return_type)
                                   
    line_indx_Hbeta = measure_line_index(wave,flux,line_info=line_info_Hbeta,
                                   num_refit = num_refit, Gauss=Gauss,
                                   return_type=return_type)
                                   
    line_indx_Halpha = measure_line_index(wave,flux,line_info=line_info_Halpha,
                                   num_refit = num_refit, Gauss=Gauss,
                                   return_type=return_type)
                                   
    return np.hstack([line_indx_Hdelta, line_indx_Hgamma, 
                      line_indx_Hbeta, line_indx_Halpha])


if __name__ == '__main__':
    print('')
    print('@Yayaha: testing ''measure_line_index'' ...')
    print('')
    print('@Yayaha: totally good spectrum:')
    fp1 = '/Users/yaoyuhan/Documents/data/Lamo_all_Mira/spec-56737-HD172545N271812V01_sp06-003.fits'
    l1=measure_Balmer_series(fp1, num_refit=(100,None), Gauss=True,return_type='dict')
    print l1
    print('')
    print('@Yayaha: totally bad spectrum:')
    fp2 = '/Users/yaoyuhan/Documents/data/Lamo_all_Mira/spec-56927-GAC062N19B1_sp10-046.fits'
    l2=measure_Balmer_series(fp2, num_refit=(100,None), Gauss=True,return_type='dict')
    print l2
    
    
    