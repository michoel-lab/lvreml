import numpy as np 
#import lvreml.py as py
from lvreml.modules.load_HLC_dataset import load_HLC_dataset
from lvreml.figcodes.Likelihoodplots import Likelihoodplots
from lvreml.figcodes.VariancePlot import VariancePlot
from lvreml.figcodes.residualplots import residualplots
from lvreml.figcodes.HeatmapPlots import HeatmapPlots
from lvreml.figcodes.SNPvsTheta import SNPvsTheta
from lvreml.figcodes.RhovsHidden import RhovsHidden

""" Set the location of data files """
geno_location = '/home/ammar/Downloads/curatedGenotype/genodata.npz'
expr_location = '/home/ammar/Downloads/curatedExpressionLiver/'
features_location = '/home/ammar/Downloads/curatedGenotype/features.txt'

fig_save_location = '/home/ammar/Downloads/LVREML_figs/'

""" Load the dataset """
C,Znall,Yn = load_HLC_dataset(geno_location,expr_location,features_location)

theta =  np.array([0,5,10,20]) # No. of known factors
rho = np.arange(5,120,10) # No. of hidden factors

""" Plot the variance explained by known factors """
VariancePlot(C,Znall,fig_save_location)

""" Plot the Likelihood values for given values of theta and rho """
""" flagC = True, indicates the use of LVREML with PC's as known confounders,
    (to test whether the known optimal latent variables (also principal components) are found)
 flagC = False, indicates the use of LVREML with given known confounders"""
Likelihoodplots(C,Znall,Yn,theta,rho,fig_save_location,flagC=True)
""" Plot the residual values """
residualplots(C,Znall,Yn,theta,rho,fig_save_location,flagC=True)

""" Plot the overlap matrix """
HeatmapPlots(C,Znall,Yn,theta[1:3],50,fig_save_location,flagC=True)

""" Plot No. of Known factors(SNPs) selected with respect to theta values """
SNPvsTheta(C,Znall,theta,fig_save_location)

"""  Plot No. of inferred Hidden factors(SNPs) selected with respect to rho values """
rho = np.arange(0.0,0.95,0.05)
RhovsHidden(C,Znall,Yn,theta,rho,fig_save_location)

