#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:43:00 2019

@author: ammar
"""
import numpy as np
import matplotlib.pyplot as plt
from lvreml.modules.initial_screen import initial_screen
#%matplotlib inline
import time

def SNPvsTheta(C,Znall,theta,location):
	tic = time.time()
	sel = np.empty(len(theta))
	for i in range(len(theta)):
	    beta2,varexpl,idx = initial_screen(C,Znall,theta[i])
	    sel[i] = len(idx)

	plt.figure(figsize=(10, 6), dpi= 80, facecolor='w', edgecolor='k')
	plt.plot(theta,sel,'-bo')
	plt.grid(which='major')
	plt.grid(which='minor')
	plt.minorticks_on()
	plt.ylabel('No. of SNPs Selected',fontsize=18)
	plt.xlabel('Theta',fontsize=18)
	plt.xticks(fontsize=15)
	plt.yticks(fontsize=15)
	plt.savefig(location+'SNPsvsTheta.png')
#	toc = time.time()
#	print('Elapsed time:%.2f secs'%(toc-tic))
#	plt.show()
