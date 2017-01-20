# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:49:38 2016

@author: yaoyuhan

This is from 'https://github.com/vr2262/azel2radec'
"""
import numpy as np
import numexpr
import math


def ltc2utc(year, month, day, hour, longitude_deg):
    """Convert Local Time to Gregorian time
    """
    if longitude_deg>0:
        hour -= math.ceil(longitude_deg/15.)
    elif longitude_deg<0:
        hour -= math.floor(longitude_deg/15.)
    assert min(hour)>=0. and max(hour)<=24.
    return year, month, day, hour
    

def jdcnv(year, month, day, hour):
    """Convert Gregorian time (UTC) to Julian Date
    
    Input:
    year  = year (scalar int)
    month = month 1-12 (scalar int)
    day   = day 1-31 (scalar int)
    hour  = fractional hour (scalar double)
    
    Output:
    julian = Julian Date (scalar double)
    
    Original IDL source available at:
    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/jdcnv.pro
    
    Revision history of IDL source:
    Converted to IDL from Don Yeomans Comet Ephemeris Generator,
    B. Pfarr, STX, 6/15/88
    Converted to IDL V5.0   W. Landsman   September 1997
       Added checks on valid month, day ranges W. Landsman July 2008
    """
    
    year = long(year)
    month = long(month)
    day = long(day)
    # account for leap years
    leap = long((month - 14) / 12.0)
    
    julian = day - 32075L + long((1461L) * (year + 4800L + leap) / 4.0) \
             + long((367L) * (month - 2 - leap * 12) / 12.0) \
             - long(3 * (long((year + 4900L + leap) / 100.0)) / 4.0) \
             + (hour / 24.0) - 0.5
    
    return julian



def ct2lst(year, month, day, hour, longitude_deg):
    """Convert Civil Time (as Julian date) to Local Sidereal Time
    
    Input:
    Local time
    
    Output:
    lst = Local Sidereal Time (scalar or vector double)
      
    The constants used in ct2lst come from Astronomical Algorithms by Jean
    Meeus, p. 84 (Eq. 11-4).
    
    Original IDL source available at:
    http://idlastro.gsfc.nasa.gov/ftp/pro/astro/ct2lst.pro
    
    Revision history of IDL source:
    Adapted from the FORTRAN program GETSD by Michael R. Greason, STX, 
              27 October 1988.
    Use IAU 1984 constants Wayne Landsman, HSTX, April 1995, results 
              differ by about 0.1 seconds  
    Longitudes measured *east* of Greenwich   W. Landsman    December 1998
    Time zone now measure positive East of Greenwich W. Landsman July 2008
    Remove debugging print statement  W. Landsman April 2009
    """
    hour = np.array(hour)
    year, month, day, hour = ltc2utc(year, month, day, hour, longitude_deg)
    julian_date = jdcnv(year, month, day, hour)
        
    c1 = 280.46061837
    c2 = 360.98564736629
    c3 = 0.000387933
    c4 = 38710000.0
    t0 = numexpr.evaluate('julian_date - 2451545.0')
    t = numexpr.evaluate('t0 / 36525')
    
    # Compute GST in seconds
    theta = numexpr.evaluate('c1 + (c2 * t0) + (t**2)'\
                             ' * (c3 - (t / c4))')
    
    # Compute LST in hours
    lst = numexpr.evaluate('theta + longitude_deg / 15.0')
    
    # Deal with LST out of bounds
    negative = numexpr.evaluate('lst < 0.0')
    negative_lst = lst[negative]
    negative_lst = numexpr.evaluate('24.0 + negative_lst % 24')
    lst[negative] = negative_lst
    lst = numexpr.evaluate('lst % 24.0') 
    
    return lst

