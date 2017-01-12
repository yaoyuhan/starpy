# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 21:49:38 2016

@author: yaoyuhan
"""
import os


def walk_dir(dirpath):
    filename_list = []
    for parent, dirnames, filenames in os.walk(dirpath):
        filename_list.extend([os.path.join(parent, filename)
                              for filename in filenames])
    n = len(filename_list)
    filename_list = filename_list[1:n+1]
    return filename_list  