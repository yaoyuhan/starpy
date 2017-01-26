#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:39:12 2017

@author: yaoyuhan
"""
import numpy as np


def MKreftable(path_disp, path_comp, path_ref):
    '''
    Generate reference table used in iraf refspec function.
    '''
    f_comp = open(path_comp,mode='r')
    f_disp = open(path_disp,mode='r')
    a = f_comp.readlines()
    b = f_disp.readlines()
    f_comp.close()
    f_disp.close()
    m=len(a)
    n=len(b)
    ind = (-1)*np.ones(n)
    for i in range(n):
        l = b[i][:-9]
        for j in range(m):
            if a[j][1:-6]>l:
                ind[i]=j
            if ind[i]!=-1:
                break
    f_ref = open(path_ref, mode='w')
    ind = np.int_(ind)
    for i in range(n):
        if ind[i]==0:
            f_ref.write(b[i][:-1]+' '+a[ind[i]])
        else:
            f_ref.write(b[i][:-1]+' '+a[ind[i]-1][:-1]+' '+a[ind[i]])
    f_ref.close()        
