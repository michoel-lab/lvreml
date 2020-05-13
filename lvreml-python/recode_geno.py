#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:44:11 2019

@author: ammar
"""

# Script for recoding genotype data to 0,1,2 & -1 
import numpy as np
location = '/home/ammar/Downloads/curatedGenotype/'
def recode_geno(location):
    f = open(location+'genotype.txt','r')
    tline = f.readline()
    numsample = len(tline.split()[1:]) # first column is feature_id
    for numline,l in enumerate(f):
        pass
    numline = numline+1
    data = np.ones((numline,numsample))*8
    ftrs = np.zeros((numline,),dtype='int64')
    indv = np.int64(tline.split()[1:])
    f.close()
    f = open(location+'genotype.txt','r')
    tline = f.readline()
    for k in range(numline):
        tline = f.readline()
        rw = np.array(tline.split())
        ftrs[k] = np.int64(rw[0])
        rwdt = rw[1:]
        isna = np.where(rwdt=='NA')[0]
        rwdt = np.delete(rwdt,isna,axis=0)
        u = np.unique(rwdt)
        if len(u)==2: 
            if len(list(set(u[1])))==2:
                data[k, np.where(rw==u[0])[0]-1] = 0
                data[k, np.where(rw==u[1])[0]-1] = 1
                data[k, np.where(rw=='NA')[0]-1] = -1
            if len(list(set(u[1])))==1:
                data[k, np.where(rw==u[0])[0]-1] = 1
                data[k, np.where(rw==u[1])[0]-1] = 0
                data[k, np.where(rw=='NA')[0]-1] = -1
            
        if len(u)==3 or len(u)==1:    
            data[k, np.where(rw==u[0])[0]-1] = 0
            data[k, np.where(rw=='NA')[0]-1] = -1 
            if len(u)==3:
                data[k, np.where(rw==u[1])[0]-1] = 1
                data[k, np.where(rw==u[2])[0]-1] = 2
        
  #  print(k)
    f.close()
    filename = location+'genodata.npz'
    np.savez(filename, indv=indv, ftrs=ftrs, data=data)
  #  return (indv, ftrs, data)
