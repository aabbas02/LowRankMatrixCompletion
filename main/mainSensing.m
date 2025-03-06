clc
close all
clear all
dir = pwd;
% For linux, replace '\' with '/'
cd ..
%addpath(genpath('.\functions'));
addpath(genpath('.\functionsMtrxSnsng'));
addpath(genpath('.\utils'));
cd(dir)    
%---------------------------------
altMin  = 1;
r = 5;
n = 600;
q = 1000;
m = 100;
numBlocks = 20;   %effectively, m_new = numBlocks
r_ = ones(1,numBlocks)*(m/numBlocks);
T = 50;
MC = 125;
% generate rank-r X*
Ustr = orth(randn(n,r));
Bstr = randn(r,q);
X  = Ustr*Bstr;
% generate q matrices A_k of size m x n, m << n
Ak_ = cell(q,1);
AkCllps_ = cell(q,1);
yk_ = cell(q,1);
ykCllps_ = cell(q,1);
ykPerm_ = cell(q,1);
M = zeros(n,q);
MCllps = zeros(n,q);
MPerm = zeros(n,q);
SDVals_UnPerm = zeros(MC,T+1);
SDVals_sLcl = zeros(MC,T+1);
SDVals_Perm = zeros(MC,T+1);
same = 1;
%if same
%    pi_map = get_permutation_r(m,r_);
%end
for mc = 1 : MC
    if same
        pi_map = get_permutation_r(m,r_);
    end
    for k = 1 : q
        Ak_{k} = randn(m,n);
        yk_{k} = Ak_{k}*X(:,k);
        M(:, k)  = Ak_{k}'*yk_{k};         
        if ~same
            pi_map = get_permutation_r(m,r_);
        end
        ykPerm_{k} = yk_{k}(pi_map);
        AkCllps_{k} = zeros(length(r_),n);
        ykCllps_{k} = zeros(length(r_),1);
        for i = 1 : length(r_)
            start = sum(r_(1:i)) - r_(i) + 1;
            stop  = sum(r_(1:i));
            AkCllps_{k}(i,:) = sum(Ak_{k}(start:stop,:));
            ykCllps_{k}(i) = sum(ykPerm_{k}(start:stop,:));
        end
        MCllps(:, k) = AkCllps_{k}'*ykCllps_{k};
    end
    [U0,~,~,] = svd(M,"econ");
    U0 = U0(:,1:r);
    %-----------------------------------
    [U0Cllps,~,~] = svd(MCllps,"econ");
    U0Cllps = U0Cllps(:,1:r);
    %---------------------------------------
    %[U0Perm,~,~] = svd(MPerm,"econ");
    %U0Perm = U0Perm(:,1:r);
    %------------------------------
    updtP = 0;
    SDVals_UnPerm(mc,:) = altGDMin_MtrxSensingPerm(Ak_, yk_,Ak_, yk_, U0,r,T,Ustr,r_,updtP,same,altMin);
    updtP = 1;
    SDVals_sLcl(mc,:) = altGDMin_MtrxSensingPerm(Ak_, ykPerm_, ...
    AkCllps_, ykCllps_, U0Cllps,r,T,Ustr,r_,updtP,same,altMin);
    %SDVals_Perm(mc,:) = altGDMin_MtrxSensing(Ak_, ykPerm_, U0Perm,r,T,Ustr);
    mc
end
%---
plotRslts(SDVals_sLcl, SDVals_Perm, SDVals_UnPerm,n,q,r,m,numBlocks,MC,same);
