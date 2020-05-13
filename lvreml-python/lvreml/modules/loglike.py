#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 13:15:50 2020

@author: ammar
"""
import numpy as np
def loglike(K,C):
    Ex,Vx = np.linalg.eigh(C)
    t = np.flip(np.argsort(Ex))
    Vx = Vx[:,t]
    D, V = np.linalg.eigh(K)
    Dinv = np.diag(1/D)
    t1 = np.sum(np.log(D))
    t2 = np.trace(np.dot(np.dot(np.dot(V,Dinv),V.T),C))
    LL = -(t1+t2)
    return LL