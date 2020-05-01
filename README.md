# lvreml

## Introduction

**lvreml** implements a restricted maximum-likelihood method for learning latent variance components in gene expression data. The method is fully described in the paper:

[MA Malik](https://www.uib.no/en/persons/Muhammad.Ammar.Malik) and [T Michoel](https://lab.michoel.info). *Restricted maximum-likelihood method for learning latent variance components in gene expression data with known and unknown confounders.*

A copy of the paper and supplementary information is also available in the [paper](./paper/) directory of this repository.

**lvreml** is available for Matlab and Python.

## Installation

### Matlab

Add the directory containing the m-files ([./m/](./m/)) to the Matlab search path: if you cloned or downloaded the repository to a directory **/some_path/lvreml/**, then issue the command

```matlab
>> addpath(/some_path/lvreml/m)
```

on the matlab command line.

### Python

## Usage

### Package content

Both the Matlab and Python versions contain the same functions:

* **data_prep**
* **initial_screen**
* **loglike**
* **lvreml**

#### data_prep

**data_prep** takes as input matrices of gene expression (Y) and known covariate data (Z) and performs basic sanity checks and data preparation tasks:

* Data matrices needs to be in the format where columns represent variables and rows represent samples. **data_prep** checks whether this is the case by checking along which dimension Y and Z have the same size. Beyond this, the package assumes that all data are in the right format, that is, the samples in Y and Z are aligned such that data from the same sample is in the same row in both matrices.
* Expression data is centred such that each sample has mean zero to remove systematic effects on the sample mean.
* The sample covariance matrix for the expression data is calculated.
* The covariate data are normalized to have unit L2 norm for each variable.

#### initial_screen

### Help


