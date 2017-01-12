# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:49:38 2016

@author: yaoyuhan
"""
import math
import time
from astropy.io import fits


def angle_distance(ra_1, dec_1, 
                    ra_2, dec_2):
    '''
    Input data can only be in the unit of degree.
    Output angle distance also in the unit of degree.
    (I'll imporve that later...)
    
    Return the angle distance of two coordinates.
    Refer to https://en.wikipedia.org/wiki/Angular_distance
    '''
    ra_1 = math.radians(ra_1)
    dec_1 = math.radians(dec_1)
    ra_2 = math.radians(ra_2)
    dec_2 = math.radians(dec_2)
    result = math.acos(math.sin(dec_1)*math.sin(dec_2) + 
                       math.cos(dec_1)*math.cos(dec_2)*math.cos(ra_1-ra_2))
    result = math.degrees(result)
    return result         


def degreeTohms_dms(RA,DEC):
    """
    Parameters
    ----------
    RA: 56.263561
    DEC: 56.318726
    
    Return
    ------
    ra: 034503.3
    dec: 561907
    """
    ra_sec = []
    ra_sec_ = []
    ra_h = math.floor(RA/15)
    ra_sec.append(str(ra_h)) # ra_h
    ra_m = math.floor((RA-ra_h*15)/0.25)
    ra_sec.append(str(ra_m)) # ra_m
    ra_sec.append(str(round((RA-ra_h*15-ra_m*0.25)*3600/15,1))) # ra_s
    for i in range(3):
        if len(ra_sec[i]) !=4:
            ra_sec_.append('0'+ra_sec[i])
        else:
            ra_sec_.append(ra_sec[i])    
    ra_str = ra_sec_[0][:2] + ra_sec_[1][:2] + ra_sec_[2]
    
    dec_sec = []
    dec_sec_ = []
    dec_d = math.floor(DEC)
    dec_sec.append(str(dec_d))
    dec_m = math.floor((DEC-dec_d)*60)
    dec_sec.append(str(dec_m))
    dec_sec.append(str(round(((DEC-dec_d)*60 - dec_m)*60)))
    for i in range(3):
        if len(dec_sec[i]) !=4:
            dec_sec_.append('0'+dec_sec[i])
        else:
            dec_sec_.append(dec_sec[i])        
    dec_str = dec_sec_[0][:2] + dec_sec_[1][:2] + dec_sec_[2][:2]    
    return ra_str, dec_str   


def save_coo_lamost(star_path,filename,savedir='/Users/yaoyuhan/Desktop/'):
    n = len(star_path)
    ra = []
    dec = []
    spid = []
    fiberid = []
    mjd = []
    lmjd = []
    for i in xrange(n):
        hl = fits.open(star_path[i])[0].header
        ra.append(hl['RA'])
        dec.append(hl['DEC'])
        spid.append(hl['SPID'])
        fiberid.append(hl['FIBERID'])
        mjd.append(hl['MJD'])
        lmjd.append(hl['LMJD'])
    c1 = fits.Column(name='RA', format='D', array = ra)
    c2 = fits.Column(name='DEC', format='D', array = dec)
    c3 = fits.Column(name='SPID', format='D', array = spid)
    c4 = fits.Column(name='FIBERID', format='D', array = fiberid)
    c5 = fits.Column(name='MJD', format='D', array = mjd)
    c6 = fits.Column(name='LMJD', format='D', array = lmjd)
    cols = fits.ColDefs([c1, c2, c3, c4, c5, c6])    
    tbhdu = fits.BinTableHDU.from_columns(cols)
    prihdr = fits.Header()
    prihdr['AUTHOR'] = 'Yao Yuhan'
    prihdr['DATE'] = time.strftime('%Y-%m-%d',time.localtime(time.time()))
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(savedir+filename+'.fits') 

      
def save_coo_simple(star_path,filename,savedir='/Users/yaoyuhan/Desktop/'):
    n = len(star_path)
    ra = []
    dec = []
    for i in xrange(n):
        hl = fits.open(star_path[i])[0].header
        ra.append(hl['RA'])
        dec.append(hl['DEC'])
    c1 = fits.Column(name='RA', format='D', array = ra)
    c2 = fits.Column(name='DEC', format='D', array = dec)
    cols = fits.ColDefs([c1, c2])    
    tbhdu = fits.BinTableHDU.from_columns(cols)
    prihdr = fits.Header()
    prihdr['AUTHOR'] = 'Yao Yuhan'
    prihdr['DATE'] = time.strftime('%Y-%m-%d',time.localtime(time.time()))
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(savedir+filename+'.fits')    
    

if __name__ == '__main__':
    print('')
    print('@Yayaha: testing ''degreeTohms_dms'' ...')
    hd = degreeTohms_dms(89.5476815,23.1037762)
    print('Return Value:')
    print hd