"""
INITIAL_SCREEN - Rapid screening of univariate variance component models
INITIAL_SCREEN screens known covariates for possible filtering before
applying the lvreml function. Screening is based on estimating the
variance explained by each covariate alone, followed by a pruning step to 
obtain a set of non-redundant, linearly independent set of covariates. 
Details and rationale are in Section S4 of the paper.

USAGE: [beta2,varexpl] = initial_screen(C,Z)
       [beta2,varexpl,idx] = initial_screen(C,Z,theta)
INPUT: C - (n x n) empirical sample covariance matrix of expression data, 
           see also data_prep
       Z - (n x d) matrix of normalized data for d covariates (known
           confounders) in n samples, see also data_prep
OUTPUT: beta2   - (d x 1) vector, variance parameter in univariate lvreml
                  model for each covariate, see Section S4 of the paper
        varexpl - (d x 1) vector, variance in C explained by each
                  covariate (= beta2/tr(C))
        idx     - indices for a subset of linearly independent covariates
                  with varexpl>=theta, only returned if number of input
                  arguments equals 3
                  
AUTHOR: Muhammad Ammar Malik
        muhammad.malik@uib.no
        https://ammarmalik93.github.io/
REFERENCE: MA Malik and T Michoel. Restricted maximum-likelihood method
for learning latent variance components in gene expression data with
known and unknown confounders. 
"""
import numpy as np
from numpy.linalg import matrix_rank

def initial_screen(C,Z,vcut):
    ns,nc = Z.shape # no. of covariates and samples
    beta2 = np.zeros((nc,),dtype=float)
    for k in range(nc):
        beta2[k] = (ns*(np.dot(np.dot(Z[:,k].T,C),Z[:,k]))/(ns-1)) - (np.trace(C)/(ns-1))
    varexpl = beta2/np.trace(C)
#    varexpl = varexpl[varexpl>0]
    # get best linearly independent set with varexpl>threshold, if necessary
    vs = np.sort(varexpl)[::-1]
    t = np.argsort(varexpl)[::-1]
    if vcut == 0:
        print("No covariates to be selected") 
        idx = np.array([], dtype=int)
    elif vcut < 1 and vcut > 0:
 #       print("Finding linearly independent covariates w.r.t given threshold %.4f"%vcut)   
        idx = t[vs>vcut]
        for k in range(1,len(idx)):
            ID = idx[0:k+1]
            rnk = matrix_rank(Z[:,ID[ID>0]])
            if rnk < (ID>0).sum():
                idx[k] = 0
        idx = idx[idx>0]
        
    else:
        print("Finding best %d linearly independent covariates"%vcut)
        idx = t
        temp = 1
        for k in range(1,len(idx)):
            ID = idx[0:k+1]
            rnk = matrix_rank(Z[:,ID[ID>0]])
            if rnk < (ID>0).sum():
                idx[k] = 0
            else:
                temp = temp+1
                if temp == vcut:
                    break
        idx = idx[idx>0]
        idx = idx[0:vcut]

    
    return beta2,varexpl,idx
