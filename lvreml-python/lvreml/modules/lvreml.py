#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 12:41:24 2019

@author: ammar
"""

import numpy as np
from .data_prep import data_prep

def lvreml(Y,Z,targetX):
    # LVREML - Restricted maximum-likelihood solution
    # for linear mixed models with known and latent
    # variance components
    
    C,Z,_ = data_prep(Y,Z)
    ns = len(C) # number of samples
    trC = np.trace(C)
    
    # Full SVD of known covariates
    if Z.size != 0:
        U,G,V = np.linalg.svd(Z)
        V = V.T
        nc = G.shape[0] # number of linearly independent covariates
        
        # Split in space spanned by covariates and orthogonal complement
        G1 = G*np.eye(nc)
        V1 = V[:,0:nc]
        U1 = U[:,0:nc]
        U2 = U[:, nc:]

        
        # Corresponding matrix blocks
        C11 = (U1.T).dot(C).dot(U1)
        C22 = (U2.T).dot(C).dot(U2)
        EC11,_ = np.linalg.eigh(C11)
        
        # Get REML estimate for X       
        if C22.size != 0:
            Ex,Vx = np.linalg.eigh(C22)
            Evx = np.sort(Ex)
            t = np.argsort(Ex)
            # minimum eigenvalue on known covariate space, needed to
            # check validitiy of the solution
            lambdamin = EC11.min()
            # variance not explained by X, at all possible number of 
            # latent variables; this is called "f(p)" in the appendix
            # of the paper
            sigma2vec = np.cumsum(Evx)/range(1,Evx.shape[0]+1)
#            sigma2vec = np.round(sigma2vec,4)
            # check if 3rd argument wants a given number of latent
            # variables or target variance explained; in both cases
            # adjust to make sure the solution will be valid
            if targetX >= 1:
                nx = np.round(targetX) # number of latent variables
                # find the smallest number of latent variables we must
                # include (largest number we can cut)
                ncut_max = (np.array(np.logical_and(sigma2vec<lambdamin,Evx>Evx[0]).nonzero()).T[-1][0])
                # cut fewer if we've asked for more latent variables
                if ncut_max < Evx.shape[0]-nx:
#                if nx < ncut_max:
                    print("Number of latent variables asked for is too small")
                    print("Infering minimum number of latent variables: %d"%(Evx.shape[0]-ncut_max-1))
                    ncut = ncut_max+1
                elif nx > Evx.shape[0]:
                    print("Number of latent variables asked for is too large")
                    print("Infering maximum number of latent variables: %d"%Evx.shape[0])
                    ncut = np.array([])
                else:
                    print("Inferring desired number of latent variables: %d"%(nx))
                    ncut = Evx.shape[0]-nx
 #               ncut = min([Evx.shape[0]-nx, ncut_max])
            elif targetX >= 0:
                varexpl = targetX # target variance explained
                # Set the target residual variance
                resvar = min([(1-varexpl)*trC/ns, lambdamin])
                #resvar = np.round(resvar,4)
                # find number of variables that need to be cut to explain
                # no less than the target variance explained
                ncut = (np.array(np.logical_and(sigma2vec<resvar,Evx>Evx[0]).nonzero()))#.T[-1][0]+1
                if ncut.size == 0:
                    ncut = np.array([])
                else:
                    ncut = ncut.T[-1][0]+1
                
            else:
                raise Exception('lvreml::lvreml::3rd argument must be integer or value between zero and one')
            
            if ncut.size == 0:
                ncut = 0
                sigma2 = 0
            else:
                # estimate sigma^2 parameter
                sigma2 = Evx[0:ncut+1].mean()
            
            # pull back selected eigenvectors of C22 to original space
            X = np.dot(U2,Vx[:,np.flip(t[ncut:])]);
            alpha2 = np.flip(Evx[ncut:]-sigma2)
            
        else:
            X = np.array([])
            alpha2 = np.array([])
            sigma2 = 0
        
        # Get estimates for the covaraicne matrix among the known
        # confounders
        V1G = (V1).dot(np.linalg.inv(G1))
        B = (V1G).dot(C11-sigma2*np.eye(nc)).dot(V1G.T)
        
        # Get estimates for the covariance matrix between the known
        # confounders and latent variables
        
        if X.size != 0:
            D = (V1G).dot(U1.T).dot(C).dot(X)
        else:
            D = np.array([])
        
        # Return covariance matrix K
        K = np.c_[Z,X].dot((np.r_[(np.c_[B,D]),(np.c_[D.T, np.diag(alpha2)])])).dot(np.r_[Z.T,X.T]) + np.eye(ns)*sigma2
        K = 0.5*(K+K.T)
        
    else:
        # no covariates provided, only estimating hidden components
        # using PCA
        print("no covariates provided, only estimating hidden components using PCA ")
        B = np.array([])
        D = np.array([])
        Ex,Vx = np.linalg.eigh(C)
        Evx = np.sort(Ex)
        t = np.argsort(Ex)
        Vx = Vx[:,t]
        sigma2vec = np.cumsum(Evx)/range(1,Evx.shape[0]+1)
           
#        nmin = np.array((Evx[1:]-sigma2vec[:len(sigma2vec)-1]>0).nonzero()).T[0][0]
        nmax = np.array((Evx[1:]-sigma2vec[:len(sigma2vec)-1]>0).nonzero()).T[-1][0]
        if targetX >= 1:
#            nx = np.round(targetX)
            nx = ns - np.round(targetX)
#            nx = Evx.shape[0]-np.round(targetX)
#            sigma2 = Evx[0:len(Evx)-nx].mean()
        elif targetX >= 0:
            varexpl = targetX        
            resvar = (1-varexpl)*trC/ns
            nx = np.array((sigma2vec<resvar).nonzero()).T[-1][0]+1
#            nx = ns-nx
        else:
            raise Exception('lvreml::lvreml::3rd argument must be integer or value between zero and one')
        
            
#        nx = max([nx, nmin])
#        sigma2 = Evx[0:len(Evx)-nx].mean()
#        X = np.fliplr(Vx[:,-nx:])
#        alpha2 = np.flip(Evx[-nx:]-sigma2)
        nx = min([nx,nmax])
#        print('NX =',nx)
        sigma2 = Evx[0:nx].mean()
        alpha2 = np.flip(Evx[-(len(Evx)-nx):]-sigma2)
        X = np.fliplr(Vx[:,-(len(Evx)-nx):])
        
        # Return covariance matrix K
        K = np.dot(np.dot(X,np.diag(alpha2)),X.T) +np.eye(ns)*sigma2
        K = 0.5*(K+K.T)
        
    return X,alpha2,B,D,sigma2,K
