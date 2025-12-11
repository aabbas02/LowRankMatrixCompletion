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
n = 500; q = 500; r = 5;
m = 100; numBlocks = 10;   %effectively, m_new = numBlocks
r_ = ones(1,numBlocks)*(m/numBlocks);
T = 200;
TAltMin = 25; %0.5*T+1; % Outer AltMin Iterations 
T_LS = 200; % Maximum GD iterations for each LS problem,usually terminates because of norm of gradient
MC = 1;
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
time_UnPerm = zeros(MC,T+1); SDVals_UnPerm = zeros(MC, T+1);
%-------------------------------------------------------
SDVals_sLcl = zeros(MC,T+1); time_sLcl = zeros(MC,T+1);
SDVals_AltMinExct = zeros(MC,TAltMin+1); time_AltMinExct=zeros(MC,TAltMin+1);
SDVals_AltMin = zeros(MC,TAltMin+1); time_AltMin=zeros(MC,TAltMin+1);
%-------------------------------------------------------
SDVals_sLclCllps = zeros(MC,T+1); time_sLclCllps = zeros(MC,T+1);
SDVals_AltMinExctCllps = zeros(MC,TAltMin+1); time_AltMinExctCllps = zeros(MC,TAltMin+1);
SDVals_AltMinCllps = zeros(MC,TAltMin+1); time_AltMinCllps =zeros(MC,TAltMin+1);
%------------------------------------------------------
eta_c = 0.3;
eta_L = 1;
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
    %--- Unpermuted
    %updtP = 0; altMin = 0; exact = 0;
    %[SDVals_UnPerm(mc,:),time_UnPerm(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, yk_,Ak_, yk_, U0,r, ...
    %    T,Ustr,r_,updtP,same,altMin,T_LS,exact,eta_c);
    % ---------------------
    % ----- COLLAPSED -----
    % ---------------------
    % AltGDMin with P - Collapsed Only
    updtP = 1; altMin = 0; exact = 0; cllpsOnly = 1;
    [SDVals_sLclCllps(mc,:), time_sLclCllps(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_,AkCllps_, ykCllps_, U0Cllps,r, ...
        T,Ustr,r_,updtP,same,altMin,T_LS,exact,eta_c,eta_L,cllpsOnly);
    %--- AltMin Kronecker LS with P - Collapsed Only
    %updtP = 1; altMin = 1; exact = 1; cllpsOnly = 1;
    %[SDVals_AltMinExctCllps(mc,:), time_AltMinExctCllps(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_, AkCllps_, ykCllps_, U0Cllps, ...
    %    r,TAltMin,Ustr,r_,updtP,same,altMin,T_LS,exact,eta_c,eta_L,cllpsOnly);    
    %--- AltMin using GD with P - Collapsed Only
    updtP = 1; altMin = 1; exact = 0; cllpsOnly = 1;
    [SDVals_AltMinCllps(mc,:), time_AltMinCllps(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_, AkCllps_, ykCllps_, U0Cllps, ...
        r,TAltMin,Ustr,r_,updtP,same,altMin,T_LS,exact,eta_c,eta_L,cllpsOnly);
    % ---------------------
    % --- UNCOLLAPSED -----
    % ---------------------
    updtP = 1; altMin = 0; exact = 0; cllpsOnly = 0;
    [SDVals_sLcl(mc,:), time_sLcl(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_,AkCllps_, ykCllps_, U0Cllps,r, ...
        T,Ustr,r_,updtP,same,altMin,T_LS,exact,eta_c,eta_L,cllpsOnly);
    %--- AltMin Kronecker LS with P - not Collapsed
    %updtP = 1; altMin = 1; exact = 1; cllpsOnly = 0;
    %[SDVals_AltMinExct(mc,:), time_AltMinExct(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_, AkCllps_, ykCllps_, U0Cllps, ...
    %    r,TAltMin,Ustr,r_,updtP,same,altMin,T_LS,exact,eta_c,eta_L,cllpsOnly);    
    %--- AltMin using GD with P - not Collapsed
    updtP = 1; altMin = 1; exact = 0; cllpsOnly = 0;
    [SDVals_AltMin(mc,:), time_AltMin(mc,:)] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_, AkCllps_, ykCllps_, U0Cllps, ...
        r,TAltMin,Ustr,r_,updtP,same,altMin,T_LS,exact,eta_c,eta_L,cllpsOnly);
    mc
end
%---
plotRslts(time_sLcl, SDVals_sLcl, ...
    time_UnPerm, SDVals_UnPerm, ...
    time_AltMinExct,SDVals_AltMinExct,...
    time_AltMin, SDVals_AltMin, ...
    time_sLclCllps,SDVals_sLclCllps,...
    time_AltMinExctCllps, SDVals_AltMinExctCllps,...
    time_AltMinCllps, SDVals_AltMinCllps, ...
    n,q,r,m,numBlocks,MC,same,T_LS,eta_c,eta_L);
  