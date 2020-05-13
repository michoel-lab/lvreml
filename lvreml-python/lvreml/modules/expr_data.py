#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:27:46 2019

@author: ammar
"""
import numpy as np
import pandas as pd
def expr_data(location):
    E = pd.read_table(location+'expression.txt')
    data = np.array(E.iloc[:,1:])
    indv = np.array(E.columns[1:]).astype('int')
    return(indv,data)
