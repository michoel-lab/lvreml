#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:35:40 2019

@author: ammar
"""
import numpy as np
from numpy.linalg import matrix_rank

def initial_screen(C,Z,vcut):
    ns,nc = Z.shape # no. of covariates and samples
    beta2 = np.zeros((nc,),dtype=float)
    for k in range(nc):
        beta2[k] = (ns*(np.dot(np.dot(Z[:,k].T,C),Z[:,k]))/(ns-1)) - (np.trace(C)/(ns-1))
    varexpl = beta2/np.trace(C)
#    varexpl = varexpl[varexpl>0]
    # get best linearly independent set with varexpl>threshold, if necessary
    vs = np.sort(varexpl)[::-1]
    t = np.argsort(varexpl)[::-1]
    if vcut == 0:
        print("No covariates to be selected") 
        idx = np.array([], dtype=int)
    elif vcut < 1 and vcut > 0:
 #       print("Finding linearly independent covariates w.r.t given threshold %.4f"%vcut)   
        idx = t[vs>vcut]
        for k in range(1,len(idx)):
            ID = idx[0:k+1]
            rnk = matrix_rank(Z[:,ID[ID>0]])
            if rnk < (ID>0).sum():
                idx[k] = 0
        idx = idx[idx>0]
        
    else:
        print("Finding best %d linearly independent covariates"%vcut)
        idx = t
        temp = 1
        for k in range(1,len(idx)):
            ID = idx[0:k+1]
            rnk = matrix_rank(Z[:,ID[ID>0]])
            if rnk < (ID>0).sum():
                idx[k] = 0
            else:
                temp = temp+1
                if temp == vcut:
                    break
        idx = idx[idx>0]
        idx = idx[0:vcut]

    
    return beta2,varexpl,idx
