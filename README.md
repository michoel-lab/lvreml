# lvreml

## Introduction

**lvreml** implements a restricted maximum-likelihood method for learning latent variance components in gene expression data. The method is fully described in the paper:

[MA Malik](https://www.uib.no/en/persons/Muhammad.Ammar.Malik) and [T Michoel](https://lab.michoel.info). *Restricted maximum-likelihood method for learning latent variance components in gene expression data with known and unknown confounders.*

A copy of the paper and supplementary information is also available on [arXiv](https://arxiv.org/abs/2005.02921) or [bioRxiv](https://doi.org/10.1101/2020.05.06.080648).

**lvreml** is available for Matlab and Python.

## Installation

### Matlab

Add the directory containing the m-files ([./m/](./m/)) to the Matlab search path: if you cloned or downloaded the repository to a directory **/some_path/lvreml/**, then issue the command

```matlab
>> addpath(/some_path/lvreml/m)
```

on the matlab command line.

### Python
1. Install the LVREML-Python package using: <code> pip install LVREML </code>
2. Run the notebook [LVREML_tutorial.ipynb](https://github.com/michoel-lab/LVREML-Python/blob/main/LVREML_tutorial.ipynb) to run the analysis on yeast data.
3. Source code for the LVREML-Python package can be accessed [here.](https://github.com/michoel-lab/LVREML-Python) 


## Usage

### Package content

Both the Matlab and Python versions contain the same functions:

* **data_prep**
* **initial_screen**
* **lvreml**
* **loglike**


#### data_prep

This function takes as input matrices of gene expression (Y) and known covariate data (Z) and performs basic sanity checks and data preparation tasks:

* Data matrices needs to be in the format where columns represent variables and rows represent samples. **data_prep** checks whether this is the case by checking along which dimension Y and Z have the same size. Beyond this, the package assumes that all data are in the right format, that is, the samples in Y and Z are aligned such that data from the same sample is in the same row in both matrices.
* Expression data is centred such that each sample has mean zero to remove systematic effects on the sample mean.
* The empirical sample covariance matrix for the expression data is calculated.
* The covariate data are normalized to have unit L2 norm for each variable.

The expression covariance matrix, and centred/normalized data are returned for use by the other package functions.

#### initial_screen

This function performs a rapid screening of univariate variance component models to calculate the variance explained by each known covariate on its own. Optionally, a linearly independent subset of covariates with variance explained greater than a user-defined threshold is computed. If the number of potential known covariates is very large (for instance, using genetic variants as known covariates for infer latent confounders in gene expression data), then it is recommended to perform this covariate filtering on principal components of the covariate matrix (genotype PCs in the example) and use selected covariate PCs as covariates in the model. This overcomes problems associated with including linearly independent yet highly correlated covariates in the model.

#### lvreml

This is the main package function. It takes as input matrices of gene expression (Y) and known covariate data (Z), and a target value for the total variance in Y explained by the known and latent variables combined. The data_prep function is called from within lvreml, and need not be called beforehand. Any preselection of covariates needs to be performed beforehand - the lvreml function will include all covariates in Z in the final model. The lvreml function performs basic linear algebra operations as explained in the supplementary information of the paper, and returns:

* **X**: a matrix of latent variable data (rows are samples, columns are variables), also called X in the paper.
* **alpha2**: a vector of variance explained by each latent variable; this is the diagonal of the diagonal matrix A in the paper.
* **B**: a matrix of covariances among the known covariates, also called B in the paper.
* **D**: a matrix of covariances between known and latent variables, also called D in the paper.
* **sigma2**: the squared residual variance.
* **K**: the restricted maximum-likelihood sample covariance matrix inferred by the method, also called K in the paper.

#### loglike

This function takes as input the lvreml sample covariance matrix K and the empirical sample covariance matrix C, and returns the (scaled) log-likelihood value of K.

### Further help

#### Matlab

Issue the command

```matlab
>> help function_name
```

on the Matlab command line for more detailed help and usage instructions for each package function.


