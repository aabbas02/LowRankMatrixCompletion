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
n = 500; q = 500; r = 3;
m = 100; numBlocks = 50;   %effectively, m_new = numBlocks
r_ = ones(1,numBlocks)*(m/numBlocks);
T = 200;
TAltMin = 25; %0.5*T+1; % Outer AltMin Iterations 
T_LS = 200; % Maximum GD iterations for each LS problem,usually terminates because of norm of gradient
MC = 125;
same = 1; % same permutation across columns
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
SDVals_sLcl_eta1 = zeros(MC,T+1); time_sLcl_eta1 = zeros(MC,T+1);
SDVals_sLcl_eta2 = zeros(MC,T+1); time_sLcl_eta2 = zeros(MC,T+1);
SDVals_sLcl_eta3 = zeros(MC,T+1); time_sLcl_eta3 = zeros(MC,T+1);
eta1 = 0.1;
eta2 = 0.5;
eta3 = 0.9;
SDVals_AltMin = zeros(MC,TAltMin+1); time_AltMin=zeros(MC,TAltMin+1);
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
    updtP = 1; altMin = 0; exact = 0;
    [SDVals_sLcl_eta1(mc,:), time_sLcl_eta1(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_,AkCllps_, ykCllps_, U0Cllps,r, ...
        T,Ustr,r_,updtP,same,altMin,T_LS,exact,eta1);
    %---------------------------------------
    updtP = 1; altMin = 0; exact = 0;
    [SDVals_sLcl_eta2(mc,:), time_sLcl_eta2(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_,AkCllps_, ykCllps_, U0Cllps,r, ...
        T,Ustr,r_,updtP,same,altMin,T_LS,exact,eta2);
    %---------------------------------------
    updtP = 1; altMin = 0; exact = 0;
    [SDVals_sLcl_eta3(mc,:), time_sLcl_eta3(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_,AkCllps_, ykCllps_, U0Cllps,r, ...
        T,Ustr,r_,updtP,same,altMin,T_LS,exact,eta3);
    %--- AltMin using GD with P
    updtP = 1; altMin = 1; exact = 0;
    [SDVals_AltMin(mc,:), time_AltMin(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_, AkCllps_, ykCllps_, U0Cllps, ...
        r,TAltMin,Ustr,r_,updtP,same,altMin,T_LS,exact,0);

    mc
end
%---
plotEtaRslts(eta1,time_sLcl_eta1, SDVals_sLcl_eta1, ...
             eta2, time_sLcl_eta2, SDVals_sLcl_eta2,...
             eta3, time_sLcl_eta3, SDVals_sLcl_eta3,...
                   time_AltMin, SDVals_AltMin,...
                   n,q,r,m,numBlocks,MC,same,T_LS,eta3)