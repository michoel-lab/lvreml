#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 18:43:00 2019

@author: ammar
"""
import numpy as np
import pandas as pd
from .data_prep import data_prep
import warnings
warnings.filterwarnings("ignore")

def load_datasets(geno_location,expr_location):
    df_expr = pd.read_table(expr_location)
    df_geno = pd.read_table(geno_location)
    # truncate row names to match genotype data
    df_expr.index = [x.split("-")[0] for x in df_expr.index]
    
    exprdata = df_expr.values
#    expr_indv = df_expr.index.values
    geno_data = df_geno.values
#    geno_indv = df_geno.index.values
    
    Zall = geno_data
    Y = exprdata
    
    C,Znall,Yn = data_prep(Y,Zall)

    return(C,Znall,Yn)
