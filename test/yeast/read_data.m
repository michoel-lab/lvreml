%% Read and save expression data
%
%%
expr = importdata('../../data/BYxRM_1000/elife-35471-expression.csv');
%%
data = expr.data;
genes = expr.textdata(1,2:end);
samples = expr.textdata(2:end,1);
%%
save ../../data/BYxRM_1000/elife-35471-expression.mat data genes samples

%% Read and save covariate data
%
%%
cov = importdata('../../data/BYxRM_1000/elife-35471-covariates.csv');
%%
data = cov.data;
covariates = cov.textdata(1,2:end);
samples = cov.textdata(2:end,1);
%%
save ../../data/BYxRM_1000/elife-35471-covariates.mat data covariates samples

%% Read and save genotype data
%
%%
geno = importdata('../../data/BYxRM_1000/elife-35471-genotypes.csv');
%%
data = geno.data;
markers = geno.textdata(1,2:end);
samples = geno.textdata(2:end,1);
%%
[chr,position]=strtok(markers,':');
[position,allele]=strtok(position,'_');
position = str2double(strrep(position,':',''))';
allele = strrep(allele,'_','');
%%
save ../../data/BYxRM_1000/elife-35471-genotypes.mat ...
    data markers samples chr position allele