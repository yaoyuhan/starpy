# -*- coding: utf-8 -*-
"""

Author
------
Bo Zhang

Email
-----
bozhang@nao.cas.cn

Created on
----------
- Tue Mar  8 15:26:00 2016    read_spectrum

Modifications
-------------
- Fri Jul 15 16:08:00 2016    migrate read_spectrum from lamost.py

Aims
----
- read various kinds of spectra

"""

import numpy as np
from astropy.io import fits
from astropy.table import Table, Column


def reconstruct_wcs_coord_from_fits_header(fp, filesource='lamost', dim=1):
    """ reconstruct wcs coordinates (e.g., wavelenght array) """
    # get keywords
    if filesource == 'apogee':
        hdr=fits.open(fp)[1].header
        crval = hdr['CRVAL%d' % dim]
        cdelt = hdr['CDELT%d' % dim]
    if filesource == 'lamost':
        hdr=fits.open(fp)[0].header
        crval = hdr['CRVAL%d' % dim]
        cdelt = hdr['CD1_%d' % dim]
    try:
        crpix = hdr['CRPIX%d' % dim]
    except(KeyError):
        crpix = 1
    naxis = hdr['NAXIS%d' % dim]

    # reconstruct wcs coordinates
    coord = np.arange(1 - crpix, naxis + 1 - crpix) * cdelt + crval
    return 10**coord


def read_spectrum(filepath, filesource='auto', visit=0):
    """read LAMOST/APOGEE spectrum
    
    visit: only meaningful when filesource=='apogee'
    
    Returns
    -------
    specdata: astropy.table.Table
        spectra as a table

    """
    # auto-identify the spectrum origination
    if filesource == 'auto':
        telescope = fits.open(filepath)[0].header['TELESCOP']
        if telescope[:3] == 'apo':
            return read_spectrum(filepath, filesource='apogee')
        if telescope == 'LAMOST':
            return read_spectrum(filepath, filesource='lamost')
        print filesource
    # LAMOST DR3 spectrum
    if filesource == 'lamost':
        wave=reconstruct_wcs_coord_from_fits_header(filepath, 
                                                    filesource='lamost', 
                                                    dim=1)
        flux = fits.open(filepath)[0].data[0]
        hl = fits.open(filepath)[0].header
        if 'Z' in hl.keys() and hl['Z'] != -9999:
            z = hl['Z']
            wave /= 1. + z
        flux = Column(name='flux', data=flux)
        wave = Column(name='wave', data=wave)
        return Table([wave, flux])
        
    # APOGEE spectrum
    if filesource == 'apogee':
        wave=reconstruct_wcs_coord_from_fits_header(filepath, 
                                                    filesource='apogee', 
                                                    dim=1)
        wave = Column(name='wave', data=wave)
        nvisit=fits.open(filepath)[0].header['NVISITS']
        if nvisit==1:
            flux = fits.open(filepath)[1].data
            flux_err = fits.open(filepath)[2].data
            flux_mask = fits.open(filepath)[3].data
        else:
            sub=visit+1
            flux = fits.open(filepath)[1].data[sub]
            flux_err = fits.open(filepath)[2].data[sub]
            flux_mask = fits.open(filepath)[3].data[sub]
        flux = Column(name='flux', data=flux) 
        flux_err = Column(name='flux_err', data=flux_err) 
        flux_mask = Column(name='flux_mask', data=flux_mask) 
        return Table([wave,flux,flux_err,flux_mask])


def _test_read_spectrum_lamost():
    fp1 = '/Users/yaoyuhan/Documents/data/Lamo_kiso_Mira/spec-55876-GAC_089N28_B1_sp14-044.fits'
    print(fp1)
    sp = read_spectrum(fp1)
    sp.pprint()
    
    
def _test_read_spectrum_apogee():
    fp1 = '/Users/yaoyuhan/Documents/data/Apo_Mira/apStar-r6-HD_51610.fits'
    print(fp1)
    sp = read_spectrum(fp1)
    sp.pprint()


if __name__ == '__main__':
    print('')
    print('@Yayaha: testing ''read_spectrum'' ...')
    _test_read_spectrum_lamost()
    _test_read_spectrum_apogee()

"""

#######################################
filesource:   'lamost'
#######################################



#######################################
filesource:   'apogee'
#######################################
https://data.sdss.org/datamodel/files/APOGEE_REDUX/APRED_VERS/APSTAR_VERS/TELESCOPE/LOCATION_ID/apStar.html

 HDU 0  : master header with target information
 HDU 1  : Flux (10^-17 ergs/s/cm^2/Ang)                           
 HDU 2  : Error (10^-17 ergs/s/cm^2/Ang)                          
 HDU 3  : Flag mask
 HDU 4+ : ...

"""