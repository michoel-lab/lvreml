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
geno = load('~/Projects/current/lmm-gws/data/Human_Liver_Cohort/genotype.mat');
disp(geno);
%% 
% Load the expression data:
expr = load('~/Projects/current/lmm-gws/data/Human_Liver_Cohort/expression.mat');
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

%% LVREML basic examples
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
[X,alpha2,B,D,sigma2,K]=lvreml(Yn,Z,0.5);
disp(loglike(K,C))
%%
% Run lvreml with a target number of 100 latent variables, and compute the
% log-likelihood:
[X,alpha2,B,D,sigma2,K]=lvreml(Yn,Z,100);
disp(loglike(K,C))

%% LVREML with PCs as known covariates
% As in Fig 2(A) of the paper, we will use principal components as known
% covariates to test whether the known optimal latent variables (also
% principal components) are found.
%%
% Get the PCs using SVD on the expression data
[U,S,V] = svd(Yn,'econ');
%%
% 
numKnown = [0,5,10,20]; nk = length(numKnown);
numLatent = 0:2:120; nl = length(numLatent);
ll = zeros(nl,nk);
lp = zeros(nl,nk);
tic;
for k = 1:nk 
   for l = 1:nl 
       [X,alpha2,B,D,sigma2,K]=lvreml(Yn,U(:,1:numKnown(k)),numLatent(l));
       ll(l,k) = loglike(K,C);
       lp(l,k) = size(X,2);
   end
end
rt=toc;
%%
% Note the runtime:
fprintf('Used %1.2f sec. for %d calls to lvreml and loglike, average %1.3f sec. per call.\n', ...
    rt, nl*nk, rt/(nl*nk));
%%
% Plot the log-likelihood values
plot(numLatent,ll,'.-','markersize',15)
set(gca,'fontsize',12,'box','off','linewidth',1.5)
grid
lgd = legend({'0','5','10','20'},'location','southeast');
lgd.Title.String = 'No. PCs as known covariates';
lgd.Title.FontSize = 14;
lgd.FontSize = 12;
xlabel('Number of latent variables','fontsize',16);
ylabel('Log-likelihood','fontsize',16)

%% LVREML with SNPs as known covariates
% As in Fig 2(C) of the paper, we will use selected SNPs as known
% covariates.
%%
%
numKnown = [0,5,10,20]; nk = length(numKnown);
numLatent = 0:2:120; nl = length(numLatent);
ll = zeros(nl,nk);
lp = zeros(nl,nk);
tic;
for k = 1:nk
   for l = 1:nl 
       [X,alpha2,B,D,sigma2,K]=lvreml(Yn,Z(:,1:numKnown(k)),numLatent(l));
       ll(l,k) = loglike(K,C);
       lp(l,k) = size(X,2);
   end
end
rt=toc;
%%
% Note the runtime:
fprintf('Used %1.2f sec. for %d calls to lvreml and loglike, average %1.3f sec. per call.\n', ...
    rt, nl*nk, rt/(nl*nk));
%%
% Plot the log-likelihood values. We need to be careful because the
% effective number of latent variables does not always equal the target
% number
plot(lp(:,1),ll(:,1),'.-','markersize',15)
hold on
for k=2:nk
    plot(lp(:,k),ll(:,k),'.-','markersize',15)
end
hold off
set(gca,'fontsize',12,'box','off','linewidth',1.5)
grid
lgd = legend({'0','5','10','20'},'location','southeast');
lgd.Title.String = 'No. SNPs as known covariates';
lgd.Title.FontSize = 14;
lgd.FontSize = 12;
xlabel('Number of latent variables','fontsize',16);
ylabel('Log-likelihood','fontsize',16)