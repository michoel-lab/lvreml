%% LVREML application on the Human Liver Cohort
% This notebook illustrates basic usage of the LVREML package on data from
% the Human Liver Cohort

%% Prerequisites
%
%%% Add lvreml to the path
% Add lvreml to the Matlab search path. Change the directory string to your
% local setup!
lvreml_dir = '~/Code/matlab/Statistics/lvreml/m/';
addpath(lvreml_dir);
%%% Download and process the data
% The Human Liver Cohort data is available on
% <https://www.synapse.org/#!Synapse:syn4499 Synapse>. You need to parse
% the files such that you can load the data matrices into Matlab. I store
% the data in several mat-files:
%% 
% Load the genotype data:
geno = load('../../data/Human_Liver_Cohort/genotype.mat');
disp(geno);
%% 
% Load the expression data:
expr = load('../../data/Human_Liver_Cohort/expression.mat');
disp(expr);

%% Process the data
% Here is sample code to filter SNPs and genes as described in the paper.
%%
% Find common samples between genotype and expression data:
[indv,isnp,iexpr] = intersect(geno.individual_idx,expr.individual_idx);
%%
% Select SNPs by call rate, MAF, HWE and autosomes:
[cr,maf,hwe] = snpselect(geno.data(:,isnp));
tf_snp = cr==1 & maf>=0.05 & hwe>1e-6 & geno.features.chrom(geno.feature_idx)~=0;
snps = find(tf_snp);
Zall = double(geno.data(tf_snp,isnp))';
%%
% Select transcripts by no missing data:
tf_gene = sum(isnan(expr.data(:,iexpr)),2)==0;
Y = expr.data(tf_gene,iexpr)';

%% LVREML application
%
%%
% Prepare the data:
[C,Znall,Yn]=data_prep(Y,Zall);
%%
% Perform an initial screen of covariates (SNPs), and keep only a linearly
% independent subset of SNPs that explain at least 18.8% of the variance in
% Y on their own
[beta2,varexpl,idx]=initial_screen(C,Znall,.188);
%%
% Show a histogram of the variance explained by each SNP:
histogram(varexpl)
%%
% Keep the selected set of SNPs as the set of known covariates in the
% model:
Z = Znall(:,idx);
%%
% Run lvreml with a target variance explained of 50%, and compute the
% log-likelihood:
[X,alpha2,B,D,sigma2,K]=lvreml(Yn,Z,100);
disp(loglike(K,C))