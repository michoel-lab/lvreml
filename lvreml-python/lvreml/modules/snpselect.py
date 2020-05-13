#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:06:23 2019

@author: ammar
"""
import numpy as np
from scipy.stats import chi2

def snpselect(data):
    num_samp = data.shape[1]
    ncall = np.sum(data>-1,axis=1, dtype='float')
    cr = ncall/num_samp #call rate
    nAA = np.sum(data==0,axis=1)
    nAB = np.sum(data==1,axis=1)
    nBB = np.sum(data==2,axis=1)
    fA = (nAA + 0.5*nAB)/ncall;
    fB = (nBB + 0.5*nAB)/ncall;
    maf = np.min(np.array([fA,fB]),axis=0) #Minor Allele Frequency
    #Chi square test for HWE
    eAA = ncall*(fA)**2;
    eAB = 2*ncall*fA*fB;
    eBB = ncall*(fB)**2;
    X = ((nAA - eAA)**2)/eAA + ((nAB - eAB)**2)/eAB + ((nBB - eBB)**2)/eBB;
    hwe = 1 - chi2.cdf(X,1)
    return cr,maf,hwe