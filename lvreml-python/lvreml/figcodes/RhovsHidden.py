#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 12:19:44 2019

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

def RhovsHidden(C,Znall,Yn,theta,rho,location):
    THET = np.empty((len(theta),len(rho)))
    HID = np.empty((len(theta),len(rho)))
    RHO = np.empty((len(theta),len(rho)))

    for i in range(len(theta)):
        beta2,varexpl,idx = initial_screen(C,Znall,theta[i])
        Z = Znall[:,idx]
    
        for j in range(len(rho)):
            X,_,_,_,_,_ = lvreml(Yn,Z,rho[j])
            HID[i,j] = X.shape[1]
            THET[i,j] = theta[i]
            RHO[i,j] = rho[j]

    df1 = pd.DataFrame(THET.flatten()) 
    df2 = pd.DataFrame(HID.flatten())
    df3 = pd.DataFrame(np.round(RHO,3).flatten())
    df = pd.concat([df1, df2,df3], axis=1)
    fig, ax = plt.subplots(figsize=(25,15))
    labels = [x + '-known' for x in theta.astype(str)]
    sns.set(style="whitegrid")
    df.columns = ['theta','hidden','rho']
    ax = sns.barplot(x='rho', y='hidden', hue='theta', data=df)
    ax.tick_params(axis='both',labelsize='20')
    ax.xaxis.label.set_fontsize(25)
    ax.yaxis.label.set_fontsize(25)
    ax.set_title('No. of Hidden Confounders vs RHO', fontsize=25)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h, labels, fontsize=20)
    plt.grid(which='major')
    plt.grid(which='minor')
    plt.minorticks_on()
    plt.savefig(location+'NumHiddenvsRho.png')
