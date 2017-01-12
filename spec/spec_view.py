# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:49:38 2016

@author: yaoyuhan
"""
import numpy as np
import matplotlib.pyplot as plt
from starpy.spec.read_spectrum import read_spectrum
from starpy.spec.spec import Spec
from starpy.tools.walk_dir import walk_dir


def _calculate_wave_offset(n_chunks, wave_intervals=None,
                           xtick_gap_fraction=0.05, plot_start=0.):
    # total span of spec chunks
    wave_totalspan = np.sum(np.diff(wave_intervals))
    # gap is a fraction of total span
    if n_chunks > 1:
        wave_gap = xtick_gap_fraction * wave_totalspan / (n_chunks - 1)
    else:
        wave_gap = 0

    # calculate wave_offset, i.e.,
    # how much wavelength should shift (to the left side)
    wave_offset_intercept = wave_intervals[0, 0] - plot_start
    wave_offset_spec = np.hstack((np.array(0), (wave_intervals[1:,0]-wave_intervals[:-1,1]).flatten()))
    wave_offset_gap = np.roll(np.ones(n_chunks) * wave_gap, 1)
    wave_offset_gap[0] = 0
    wave_offset = np.cumsum(wave_offset_spec - wave_offset_gap) + wave_offset_intercept
    return wave_offset


def _generate_chunk_xtick_pos_lab(n_xtick_l,
                                  n_xtick_r,
                                  wave_center,
                                  xtick_pos_step=50.,
                                  xtick_lab_step=50.,
                                  wave_offset_chunk=0.,
                                  xtick_format_str=('%.0f', '%.1f')):
    # relative xtick position
    xtick_pos_rel = np.arange(-n_xtick_l, n_xtick_r+1) * xtick_pos_step
    xtick_pos = xtick_pos_rel + wave_center - wave_offset_chunk
    # xtick labels
    xtick_lab_rel = np.arange(-n_xtick_l, n_xtick_r+1) * xtick_lab_step
    xtick_lab = [xtick_format_str[0] % xtick_lab_rel_ for xtick_lab_rel_ in xtick_lab_rel]
    xtick_lab[n_xtick_l] = xtick_format_str[1] % wave_center
    return xtick_pos, xtick_lab


def _xtick_pos_lab(n_chunks,
                   wave_centers, wave_intervals, wave_offset,
                   xtick_step,
                   xtick_format_str=('%.0f', '%.1f')):
    # reshape wave_centers if necessary
    wave_centers = np.array(wave_centers)
    if wave_centers.ndim == 1:
        wave_centers = wave_centers.reshape((wave_centers.shape[0], 1))

    # initialize results
    xtick_pos, xtick_lab = [], []

    # array
    xtick_pos_step = np.ones_like(wave_centers) * xtick_step
    # float
    xtick_lab_step = xtick_step

    for i_chunk in xrange(n_chunks):
        n_xtick_l = np.int(np.abs((wave_centers[i_chunk] - wave_intervals[i_chunk][0]) / xtick_pos_step[i_chunk]))
        n_xtick_r = np.int(np.abs((wave_centers[i_chunk] - wave_intervals[i_chunk][1]) / xtick_pos_step[i_chunk]))
        xtick_pos_, xtick_lab_ = _generate_chunk_xtick_pos_lab(
            n_xtick_l,
            n_xtick_r,
            wave_centers[i_chunk],
            xtick_pos_step=xtick_pos_step[i_chunk],
            xtick_lab_step=xtick_lab_step,
            wave_offset_chunk=wave_offset[i_chunk],
            xtick_format_str=xtick_format_str)
        xtick_pos.extend(xtick_pos_)
        xtick_lab.extend(xtick_lab_)
    return xtick_pos, xtick_lab

    
    
def multi_spec_view(fp_list, wave_interval,
                    wave_lines, norm_method='max',
                    offset_perspec = 1.0, out_format='wavelength',
                    spec_list=None):
    '''
    Imput data (wave_interval, wave_lines) in the scale of wavelength.
    
    wave_lines: x*y nparray, x<=8
    '''
    # Default colors
    color_lines = ['gold', 'b', 'g', 'darkviolet',
                   'c', 'darkorange', 'hotpink', 'k']    
    color_flux = ['r','white']
    
    # Check
    n_lines = len(wave_lines)
    assert n_lines<=8
    
    # Preparation
    wave_interval = np.array(wave_interval)
    wave_lines = np.array(wave_lines)  
    n_each = []
    for i in xrange(n_lines):
        n_each.append(len(wave_lines[i]))
    n_each = np.array(n_each)
    
    # Shift wavelength_format if necessary
    # Do it after spec_chuck
    if out_format=='wavenumber':
        wave_interval = [[1./wave_interval[0][1]*10**8, 
                          1./wave_interval[0][0]*10**8]]
        for i in range(n_lines):
            for j in range(n_each[i]):
                wave_lines[i][j] = 1./wave_lines[i][j]*10**8   
    
    # Extrack spec data
    fp_list = np.array(fp_list)
    n = len(fp_list)  
    if spec_list == None:
        spec_list = []         
        for i in xrange(n):
            spec = Spec(read_spectrum(fp_list[i]))
            spec_list.append(spec) 
    if out_format == 'wavenumber':
        for i in xrange(n):
            spec_list[i].shift_to_wavenumber()
    spec_chunk_list = [spec.extract_chunk_wave_interval(wave_interval)
                       for spec in spec_list]                         
    for i in xrange(n):
        if norm_method=='max':
            spec_chunk_list[i][0].norm_spec_max()
        elif norm_method=='median':
            spec_chunk_list[i][0].norm_spec_median()
    
    # Before spectra
    offset_specs = np.arange(n) * offset_perspec
    fig = plt.figure('spec_view', figsize=(20, 30))
    ax = fig.add_subplot(111)
    
    for i in range(n):
        spec_chunk = spec_chunk_list[i][0]
        ax.plot(spec_chunk['wave'],
                spec_chunk['flux'] + offset_specs[i],color_flux[0])
        if 'flux_mask' in spec_chunk.keys():
            flux_mask = spec_chunk['flux_mask']
            ind = np.all([flux_mask>1],axis=0)
            ax.plot(spec_chunk['wave'][ind],
                    spec_chunk['flux'][ind] + offset_specs[i],
                    color = color_flux[1],alpha=0.7)
        
    # Plot lines
    bottom = min(spec_chunk_list[0][0]['flux'])
    up = max(spec_chunk_list[n-1][0]['flux'])+ offset_specs[n-1]
    for i in xrange(n_lines):
        for j in xrange(n_each[i]):
            wv = wave_lines[i][j]
            ax.vlines(wv,bottom,up,color=color_lines[i])
    
    if out_format=='wavelength':
        ax.set_xlabel(r'$\lambda (\AA)$')
    else:
        ax.set_xlabel(r'$ Wavenumber (cm^{-1})$')
    ax.set_ylabel('$Flux$')
    
    
if __name__ == '__main__':
    dir_name = '/Users/yaoyuhan/Documents/data/Apo_Mira'
    fp_list = walk_dir(dir_name)[30:]
    wave_interval = [[15100., 17000.]]
    wave_lines = [[16723.524, 16755.140, 16767.939],
                  [15581.633, 15779.465, 15981.867, 16188.925,
                   16400.808, 16617.707, 16839.740],
                  [15745.017, 15753.291, 15770.150],
                  [15964.424, 16685.327],
                  [16883.701]]
    print ''
    print('@Yayaha: start to test multi_spec_view ...')
    multi_spec_view(fp_list[:10], wave_interval, 
                    wave_lines, norm_method='median',
                    offset_perspec = 2.0, out_format='wavenumber')