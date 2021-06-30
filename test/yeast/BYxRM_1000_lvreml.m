%% LVREML application on data from 1,012 yeast segregants (Albert et al. 2020)
% This notebook illustrates basic usage of the LVREML package on data from
% 1,012 yeast segregants from Albert et al, eLife 2020.

%% Prerequisites
%
%%% Add lvreml to the path
% Add lvreml to the Matlab search path. Change the directory string to your
% local setup!
lvreml_dir = '~/Code/matlab/Statistics/lvreml/m/';
addpath(lvreml_dir);

%% Load data
%
expr = load('../../data/BYxRM_1000/elife-35471-expression.mat');
geno = load('../../data/BYxRM_1000/elife-35471-genotypes.mat');
covrs = load('../../data/BYxRM_1000/elife-35471-covariates.mat');
% Set dimensions
ng = length(expr.genes);
ns = length(expr.samples);

%% Regress cell measurements out of expression via linear model
% The data above already contains a corrected expression matrix Y. This
% code shows how Y was obtained.
Y = zeros(size(expr.data));
for k=1:ng
    mdl = fitlm(covrs.data,expr.data(:,k),'CategoricalVars',[1]);
    Y(:,k) = mdl.Residuals.Raw;
end

%% Get principal components of genotype data
% These will be used as genotype/population structure covariates.
% Note that PCA centers variables by default
% Note that the covariates (Zn) are automatically scaled to have unit L2 norm
[Zn,S,V] = svd(geno.data-mean(geno.data),'econ');
vZ = diag(S.^2)./sum(diag(S.^2));

%% Center expression data and get sample covariance matrix
%
[C,Yn] = data_prep(expr.Y,[]);

%% Compute variance explained by genotype PCs in univariate models
% We don't need option to compute linearly independent subsets because all
% subsets of PCs are linearly independent by definition
[beta2,varexpl]=initial_screen(C,Zn);
[~,vx] = sort(varexpl,'descend');

%% Generate Figure 1A
% This figure shows the variance explained by the genotype PCs in the
% expression data vs their variance explained in the genotype data
figure(1)
scatter(vZ,varexpl,75,'markeredgecolor','k',...
    'markerfacecolor',[0 0.4470 0.7410])
set(gca,'linewidth',1.5,'fontsize',14,'box','off')
xlabel('Genotype variance explained','fontsize',16)
ylabel('Gene expression variance explained','fontsize',16)

%% Generate Figure 1B
%
numcov = 0:5:40;
rho = 0:0.05:0.95;
numX = zeros(length(rho),length(numcov));
realVX = zeros(length(rho),length(numcov));
for k=1:length(numcov)
    Zk = [];
    if k>1
       Zk = Zn(:,vx(1:numcov(k)));
    end
    for m=1:length(rho)
        disp([k,m]);
        [X,alpha2,B,D,sigma2,K]=lvreml(C,Zk,rho(m));
        numX(m,k) = size(X,2);
        realVX(m,k) = (sum(diag(Zk*B*Zk'))+sum(diag(X*diag(alpha2)*X')))/sum(diag(C));
    end
end
%% 
% Fig 1B
figure(2)
nn = 13;
bar(rho(3:nn),numX(3:nn,1:5));
lgd = cell(size(numcov));
for k=1:length(numcov)
   lgd{k} = num2str(numcov(k)); 
end
ll = legend(lgd,'location','west','fontsize',14);
ll.Position=[0.15 0.63 0.0608 0.2078];
text(0.09,38,'Nunber of genotype PC covariates:','fontsize',16)
set(gca,'linewidth',1.5,'fontsize',14,'box','off',...
    'xtick',rho(3:nn))
xlabel('Targeted gene expression variance explained (\rho)','fontsize',16)
ylabel('Number of inferred hidden covariates','fontsize',16)

%%
% Fig S1A
figure(2)
nn = 13;
bar(rho(1:nn),numX(1:nn,:));
lgd = cell(size(numcov));
for k=1:length(numcov)
   lgd{k} = num2str(numcov(k)); 
end
ll = legend(lgd,'location','west','fontsize',14);
ll.Position=[0.15 0.63 0.0608 0.2078];
text(0.0,675,'Nunber of genotype PC covariates:','fontsize',16)
set(gca,'linewidth',1.5,'fontsize',14,'box','off',...
    'xtick',rho(1:nn))
xlabel('Targeted gene expression variance explained (\rho)','fontsize',16)
ylabel('Number of inferred hidden covariates','fontsize',16)

%% OLD code below
%

%% Run lvreml
%
[X,alpha2,B,D,sigma2,K]=lvreml(Yn,Z,10);
%%
K2=Z*B*Z'+X*diag(alpha2)*X'+sigma2*eye(size(Yn,1));
K2=0.5*(K2+K2');
%%
%
Z = [];
[X,alpha2,B,D,sigma2,K]=lvreml(Yn,Z,.5);
%% variances explained
vx = [sum(diag(Z*B*Z')), sum(diag(X*diag(alpha2)*X')), ns*sigma2]/sum(diag(C))

%% plot varexpl for all SNPs
[u,m,n]=unique(geno.chr);
u = u([1:4,6:9,5,10:end]);
tf1 = ismember(n,6:9);
tf2 = n==5;
n(tf1) = n(tf1)-1;
n(tf2) = 9;
a = accumarray(n,geno.position,[length(u),1],@max);
chr_start = ones(size(a));
chr_start(2:end) = a(2:end);
chr_start = cumsum(chr_start);
chr_pos = chr_start(n) + geno.position;