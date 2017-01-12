# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:49:38 2016

@author: yaoyuhan
"""
import numpy as np
from starpy.spec.read_spectrum import read_spectrum
from starpy.spec.line_indice import measure_line_index
from starpy.tools.walk_dir import walk_dir
from starpy.spec.template_fit import best_template_subtraction, MCgiant_tp_list_init
  
    
def measure_C2_TiO(filepath):
    line_info_C24737 = {'line_center':         4737,
                        'line_range':          (4625, 4747),
                        'line_shoulder_left':  (4620, 4625),
                        'line_shoulder_right': (4750, 4775)}
    line_info_C25635 = {'line_center':         5400,
                        'line_range':          (5350, 5640),
                        'line_shoulder_left':  (5300, 5350),
                        'line_shoulder_right': (5650, 5700)}
    line_info_TiO6191 = {'line_center':         6191,
                         'line_range':          (6191.375, 6273.875),
                         'line_shoulder_left':  (6068.375, 6143.375),
                         'line_shoulder_right': (6374.375, 6416.875)}
    num_refit=(None,None)             
    spec = read_spectrum(filepath)  
    wave = spec['wave']
    flux = spec['flux']      
    line_indx_C24737 = measure_line_index(wave,flux,line_info=line_info_C24737,
                                   num_refit = num_refit, Gauss=False)
                                   
    line_indx_C25635 = measure_line_index(wave,flux,line_info=line_info_C25635,
                                   num_refit = num_refit, Gauss=False)
                                   
    line_indx_TiO6191 = measure_line_index(wave,flux,line_info=line_info_TiO6191,
                                   num_refit = num_refit, Gauss=False)
    return np.hstack([line_indx_C24737, line_indx_C25635, line_indx_TiO6191])
    
    
def measure_Balmer_series(filepath, 
                          H3_info=None, 
                          H4_info=None, 
                          num_refit=(None,None),
                          Gauss=False,
                          return_type='dict'):
    '''
    single path
    '''
    line_info_H6 = {'line_center':         4102,
                    'line_range':          (4096.00, 4109.00),
                    'line_shoulder_left':  (4080.00, 4094.00),
                    'line_shoulder_right': (4112.00, 4127.00)}
                        
    line_info_H5 = {'line_center':         4340,
                    'line_range':          (4337.00, 4347.00),
                    'line_shoulder_left':  (4320.00, 4335.00),
                    'line_shoulder_right': (4348.57, 4351.88)}
                        
    line_info_H4 = {'line_center':         4861,
                    'line_range':          (4857.00, 4868.00),
                    'line_shoulder_left':  (4845.00, 4855.00),
                    'line_shoulder_right': (4875.00, 4880.00)}                        
                        
    line_info_H3 = {'line_center':         6563,
                    'line_range':          (6548.00, 6578.00),
                    'line_shoulder_left':  (6505.00, 6540.00),
                    'line_shoulder_right': (6580.00, 6610.00)}                    
           
    spec = read_spectrum(filepath)  
    wave = spec['wave']
    flux = spec['flux']      
    line_indx_H6 = measure_line_index(wave,flux,line_info=line_info_H6,
                                      num_refit = num_refit, Gauss=Gauss,
                                      return_type=return_type)
                                   
    line_indx_H5 = measure_line_index(wave,flux,line_info=line_info_H5,
                                      num_refit = num_refit, Gauss=Gauss,
                                      return_type=return_type)
    
    if H4_info != None:
        line_indx_H4 = measure_line_index(wave=H4_info['wave_region'],
                                          flux=H4_info['flux_est'],
                                          line_info=line_info_H4,
                                          num_refit = num_refit, Gauss=Gauss,
                                          return_type=return_type)
    else:
        line_indx_H4 = measure_line_index(wave,flux,line_info=line_info_H4,
                                          num_refit = num_refit, Gauss=Gauss,
                                          return_type=return_type)
    
    if H3_info != None:
        line_indx_H3 = measure_line_index(wave=H3_info['wave_region'],
                                          flux=H3_info['flux_est'],
                                          line_info=line_info_H3,
                                          num_refit = num_refit, Gauss=Gauss,
                                          return_type=return_type)                               
    else:
        line_indx_H3 = measure_line_index(wave,flux,line_info=line_info_H3,
                                          num_refit = num_refit, Gauss=Gauss,
                                          return_type=return_type)                                   
    return np.hstack([line_indx_H6, line_indx_H5, 
                      line_indx_H4, line_indx_H3])
                      

def loopfun_Balmer(fp_list, Halpha_info=None, Hbeta_info=None,
                   num_refit=(None,None), Gauss=False, return_type='array'):
    '''
    lots of path
    '''  
    n=len(fp_list)
    
    data = []
    for i in range(n):
        if Halpha_info != None:
            H3_info = Halpha_info[i]
        else:
            H3_info = None
        if Hbeta_info != None:
            H4_info = Hbeta_info[i]
        else:
            H4_info = None
        d = measure_Balmer_series(fp_list[i], H3_info=H3_info,
                                  H4_info=H4_info, num_refit=num_refit,
                                  Gauss=Gauss, return_type=return_type)
        data.append(d)
    return data


def _test_loopfun_Balmer():
    fp_list = walk_dir('/Users/yaoyuhan/Documents/data/Lamo_kiso_Mira/')[0:3]
    n=len(fp_list)
    wv_tp_list, fl_tp_list = MCgiant_tp_list_init()
    region_info_Halpha = {'region_center':  6563.,
                          'region_left':    (6400., 6555.),
                          'region_right':   (6575., 6700.)}   
    Halpha_info = []
    for i in range(n):
        wv_tg = read_spectrum(fp_list[i])['wave']
        fl_tg = read_spectrum(fp_list[i])['flux'] 
        info = best_template_subtraction(wv_tg, fl_tg, wv_tp_list, fl_tp_list, 
                                         region_info=region_info_Halpha,
                                         star_type = 1)   
        Halpha_info.append(info)
    data = loopfun_Balmer(fp_list, Halpha_info=Halpha_info, 
                          Hbeta_info=None,
                          num_refit=(None,None), Gauss=False, 
                          return_type='array')
    return data
    

if __name__ == '__main__':
    print('')
    print('@Yayaha: testing ''masure_Balmer_series'' ...')
    fp='/Users/yaoyuhan/Documents/data/Lamo_kiso_Mira/spec-57457-KP061450N190923V01_sp03-138.fits'
    result=measure_Balmer_series(fp,
                                 num_refit=(100,None),
                                 Gauss=True,
                                 return_type='array')
    print result
    print('')
    print('@Yayaha: testing ''loopfun_Balmer'' ...')
    data = _test_loopfun_Balmer()  
    print data