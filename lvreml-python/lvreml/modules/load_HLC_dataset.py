#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:43:00 2019

@author: ammar
"""
import numpy as np
import pandas as pd
#from Codes.tests.recode_geno import recode_geno
from .expr_data import expr_data
from .data_prep import data_prep
from .snpselect import snpselect
import warnings
warnings.filterwarnings("ignore")

def load_HLC_dataset(geno_location,expr_location,features_location):
    feat = pd.read_table(features_location)
    feat.loc[np.where((feat.chrom=='X') | (feat.chrom=='A'))[0],['chrom']] = -1
    feat_id = np.array(pd.to_numeric(feat['feature_id'],errors='coerce')).astype('int')
    
    GENO = np.load(geno_location)
    geno_feat = GENO['geno_feat']
    geno_indv = GENO['geno_indv']
    geno_data = GENO['geno_data']
#    geno_data = np.float16(geno_data)
#    (geno_indv, geno_feat, geno_data) = recode_geno(geno_location)
    _,geno_idx,feat_idx =  np.intersect1d(geno_feat,feat_id,return_indices=True)
    
    chrom = np.array(feat.iloc[feat_idx,1]).astype(int)
    geno_data = geno_data[geno_idx,:]
    (expr_indv, exprdata) = expr_data(expr_location)
    
    indv,isnp,iexpr = np.intersect1d(geno_indv,expr_indv,return_indices=True)
    
    cr,maf,hwe = snpselect(geno_data[:,isnp])
    
    tf_snp = (cr==1)&(maf>=0.05)&(hwe>1e-6)&(chrom!=0)
    Zall = geno_data[tf_snp,:]
    Zall = Zall[:,isnp].T
    tf_gene = np.sum(np.isnan(exprdata[:,iexpr]),axis=1)==0
    Y = exprdata[tf_gene,:]
    Y = Y[:,iexpr].T
    
    C,Znall,Yn = data_prep(Y,Zall)

    return(C,Znall,Yn)
