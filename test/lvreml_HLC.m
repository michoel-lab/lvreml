%% Individuals
%
load ../../data/Human_Liver_Cohort/individuals.mat;

%% Genotype data 
%
geno = load('../../data/Human_Liver_Cohort/genotype.mat');

%% Expression data
%
expr = load('../../data/Human_Liver_Cohort/expression.mat');

%% Select common samples, SNPs and transcripts
% 
%%
% Common samples
[indv,isnp,iexpr] = intersect(geno.individual_idx,expr.individual_idx);
%%
% Select SNPs by call rate, MAF, HWE and autosomes
[cr,maf,hwe] = snpselect(geno.data(:,isnp));
tf_snp = cr==1 & maf>=0.05 & hwe>1e-6 & geno.features.chrom(geno.feature_idx)~=0;
snps = find(tf_snp);
Zall = double(geno.data(tf_snp,isnp))';
%%
% Select transcripts by no missing data
tf_gene = sum(isnan(expr.data(:,iexpr)),2)==0;
Y = expr.data(tf_gene,iexpr)';

%%
%
[C,Znall,Yn]=data_prep(Y,Zall);
[beta2,varexpl,idx]=initial_screen(C,Znall,.19);
Z = Znall(:,idx);
[X,alpha2,B,sigma2]=lvreml(Yn,Z,0.5);
K=Z*B*Z'+X*diag(alpha2)*X'+sigma2;

%% SNP-Gene associations
%
%[I,J,P] = kruX(Y',Zall',1e-3,1e3);
% %%
% % Save
% save ../../results/Human_Liver_Cohort/HLC_kruX_lvreml.mat I J P
%%
% Make sparse matrix
Pmat = sparse(I,J,-log10(P),size(Y,2),size(Zall,2));
%%
% Number of significant genes per SNP
nsig = sum(Pmat~=0)';
ssig = sum(Pmat)';