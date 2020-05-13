#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:43:00 2019

@author: ammar
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def VariancePlot(C,Znall,location):
	fig, ax = plt.subplots(figsize=(10,6))
	ns,nc = Znall.shape # no. of covariates and samples
	beta2 = np.zeros((nc,),dtype=float)
	for k in range(nc):
	    beta2[k] = (ns*(np.dot(np.dot(Znall[:,k].T,C),Znall[:,k]))/(ns-1)) - (np.trace(C)/(ns-1))
	varexpl = beta2/np.trace(C)
	plt.grid(True)
	plt.xlabel('Variance Explained',fontsize=18)
	plt.xticks(fontsize=15)
	plt.yticks(fontsize=15)
	sns.distplot(varexpl)
	#plt.show()
	plt.savefig(location+'VarianceExplainedBySNPs.png')

