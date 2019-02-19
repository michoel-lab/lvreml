function [callrate,maf,hwe] = snpselect(data)
% SNPSELECT - Compute call rate, MAF and HWE p-value for a SNP data set
% We assume Birdseed output file format, i.e.
%   -1 = missing call
%    0 = AA
%    1 = AB
%    2 = BB
%
%
% Copyright 2011, Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de


% number of individuals
np = size(data,2);

% call rate
ncall = sum(data>-1,2);
callrate = ncall/np;

AA = (data==0);
AB = (data==1);
BB = (data==2);

nAA = sum(AA,2);
nAB = sum(AB,2);
nBB = sum(BB,2);

fA = (nAA + 0.5*nAB)./ncall;
fB = (nBB + 0.5*nAB)./ncall;

% minor allele frequency
maf = min(fA,fB);

% chi-square test for HWE
eAA = ncall.*(fA).^2;
eAB = 2*ncall.*fA.*fB;
eBB = ncall.*(fB).^2;

X = (nAA - eAA).^2./eAA + (nAB - eAB).^2./eAB + (nBB - eBB).^2./eBB;
hwe = 1-chi2cdf(X,1);