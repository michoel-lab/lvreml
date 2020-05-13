#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 13:03:12 2019

@author: ammar
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from lvreml.modules.initial_screen import initial_screen
from lvreml.modules.lvreml import lvreml
import warnings
warnings.filterwarnings("ignore")

def residualplots(C,Znall,Yn,theta,rho,location,flagC):
    RES = np.empty((len(theta),len(rho)))
    HID = np.empty((len(theta),len(rho)))

    plt.figure(figsize=(10, 6), dpi= 80, facecolor='w', edgecolor='k')
    col = ['-bo','--rx','-.ks',':go']
    Ex,Vx = np.linalg.eigh(C)
    t = np.flip(np.argsort(Ex))
    Vx = Vx[:,t]
    for i in range(len(theta)):
        beta2,varexpl,idx = initial_screen(C,Znall,theta[i])
        if flagC == False:
            Z = Znall[:,idx]
        else:
            Z = Vx[:,0:theta[i]]
        for j in range(len(rho)):
            X,alpha2,B,D,sigma2,K = lvreml(Yn,Z,rho[j])
            RES[i,j] = sigma2
            HID[i,j] = X.shape[1]
        plt.plot(HID[i,:],RES[i,:],col[i])
    
 #   print(RES)
    labels = [x + '-known' for x in theta.astype(str)]  
    plt.grid(which='major')
    plt.grid(which='minor')
    plt.minorticks_on()
    plt.ylabel('Residual Values',fontsize=18)
    plt.xlabel('No. of Hidden confounders',fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(labels, fontsize=15)
    if flagC == False:
        plt.savefig(location+'ResvsNumHidd_LVERML.png')
    else:
        plt.savefig(location+'ResvsNumHidd_LVERML(usingC).png')
#    plt.show()
