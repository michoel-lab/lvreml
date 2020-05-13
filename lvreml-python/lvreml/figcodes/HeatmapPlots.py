#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:43:00 2019

@author: ammar
"""
import numpy as np
import matplotlib.pyplot as plt
from lvreml.modules.initial_screen import initial_screen
from lvreml.modules.lvreml import lvreml
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

def HeatmapPlots(C,Znall,Yn,theta,rho,location, flagC):
    fig, axs = plt.subplots(1,2, figsize=(20,20))
    sub = np.array([axs[0],axs[1]])
    labels = [x + '-known' for x in theta.astype(str)]
    Ex,Vx = np.linalg.eigh(C)
    t = np.flip(np.argsort(Ex))
    Vx = Vx[:,t]
    for i in range(len(theta)):
        beta2,varexpl,idx = initial_screen(C,Znall,theta[i])
        if flagC == False:
            Z = Znall[:,idx]
        else:
            Z = Vx[:,0:theta[i]]
        X,_,_,_,_,_ = lvreml(Yn,Z,rho)
#        Corr = np.dot(X.T,X)
        if theta[i]==0:
            Corr = np.dot(X.T,X)
        else:
            Corr = np.dot(X.T,Z)
            
        sub[i].set_title(labels[i], fontsize=10)
        sns.heatmap(Corr, cmap='bwr',vmin=-1, vmax=1, center=0,square= True, ax=sub[i])
#    plt.show()
    if flagC == False:
        plt.savefig(location+'CorrelationMatrix_LVERML_HiddenvsKnown.png')
    else:
        plt.savefig(location+'CorrelationMatrix_LVERML_HiddenvsKnown(usingC).png')
