#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:14:26 2019

@author: ammar
"""
import numpy as np
import warnings
def data_prep(Y,Z):
    if Z.size != 0:
        if Y.shape[0] != Z.shape[0]:
            # rows are not matching, check columns
            if Y.shape[1] == Z.shape[1]:
                # columns are samples -> transpose both
                Y = Y.T
                Z = Z.T
            else:
                raise Exception('Cannot determine sample dimension')
    else:
        warnings.warn("No covariates provided, assuming expression data is in format samples x genes")
    
    # Center Y to remove fixed effects on mean and get overlap matrix         
    Yn = Y - np.array(Y.mean(axis=1)).reshape(Y.shape[0],1)
    C = np.dot(Yn,Yn.T)/Y.shape[1]
    
    if Z.size != 0:
        # Normalize covariate to have unit L2-norm
#        Zn = (Z.T/np.sqrt(np.sum(Z**2,axis=0)).reshape(Z.shape[1],1)).T
        Zn = (Z.T/np.sqrt(np.sum(np.square(Z),axis=0)).reshape(Z.shape[1],1)).T
    else:
        Zn = np.empty((0,0))
    return C,Zn,Yn
