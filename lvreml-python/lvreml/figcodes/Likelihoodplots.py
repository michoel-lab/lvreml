#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:43:00 2019

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

def Likelihoodplots(C,Znall,Yn,theta,rho,location, flagC):
    THET = np.empty((len(theta),len(rho)))
    LL = np.empty((len(theta),len(rho)))
    HID = np.empty((len(theta),len(rho)))
    RHO = np.empty((len(theta),len(rho)))
    Ex,Vx = np.linalg.eigh(C)
    t = np.flip(np.argsort(Ex))
    Vx = Vx[:,t]
    plt.figure(figsize=(10, 6), dpi= 80, facecolor='w', edgecolor='k')
    col = ['-bo','--rx','-.ks',':go']
    
    for i in range(len(theta)):
        beta2,varexpl,idx = initial_screen(C,Znall,theta[i])
        if flagC == False:
            Z = Znall[:,idx]
        else:
            Z = Vx[:,0:theta[i]]

        for j in range(len(rho)):
            X,alpha2,B,D,sigma2,K = lvreml(Yn,Z,rho[j])   
            D, V = np.linalg.eigh(K)
            Dinv = np.diag(1/D)
            t1 = np.sum(np.log(D))
            t2 = np.trace(np.dot(np.dot(np.dot(V,Dinv),V.T),C))
            LL[i,j] = -(t1+t2)
            HID[i,j] = X.shape[1]
            THET[i,j] = theta[i]
            RHO[i,j] = rho[j]
        plt.plot(HID[i,:],LL[i,:],col[i])
	    
    labels = [x + '-known' for x in theta.astype(str)]    
    	#plt.grid(which='major')
    plt.grid()
    #	plt.minorticks_on()
    plt.ylabel('Log-likelihood value',fontsize=18)
    plt.xlabel('No. of Hidden confounders',fontsize=18)
    plt.xticks(rho,fontsize=10)
    plt.yticks(fontsize=15)
    plt.legend(labels, fontsize=15)
    if flagC == False:
    #    plt.savefig('/home/ammar/lvreml-paper-2018/notebooks/figures/LVREML/LLvsNumHidd_LVERML.png')
        plt.savefig(location+'LLvsNumHidd_LVERML.png')
    else:
        plt.savefig(location+'LLvsNumHidd_LVERML(usingC).png')
#    plt.show()

