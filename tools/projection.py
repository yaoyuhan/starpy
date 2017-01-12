# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 09:11:32 2016

@author: yaoyuhan
"""
import math
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord


def aitoff_galactic(RA,DEC,marker_size = 3.5):
    '''
    Input coordinates must be in the unit of degree.
    '''
    n = len(RA)
    c_icrs = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='icrs')
    c_galactic = c_icrs.galactic  
    l = c_galactic.l.rad
    for i in range(n):
        if l[i]>math.pi:
            l[i]=-(2*math.pi-l[i])
    b = c_galactic.b.rad
    fig = plt.figure('spec_view', figsize=(16, 10))
    ax = fig.add_subplot(111, projection = "aitoff")
    ax.set_title('Galactic projection')
    ax.grid(True)
    ax.plot(l,b, 'ro',markersize=marker_size)