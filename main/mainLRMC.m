clc
close all
clear all
dir = pwd;
% For linux, replace '\' with '/'
cd ..
addpath(genpath('.\functionsMtrxCmpltn'));
addpath(genpath('.\utils'));
cd(dir)    
tic
%---------------------------------
r = 5;
real = 0;
T_init = 5000;
T_AltMin = 1000;
n = 1000;
q = 1000;
p = 0.2;
same = 0;
MC = 15;
%------------------------
if real
    [Ustr,X,p] = getMovieLens(r);
    n = size(X,1);
    q = size(X,2);
    MC = 1;
else
    % generate rank-r X*
    Ustr = orth(randn(n,r));
    Bstr = randn(r,q);
    X  = Ustr*Bstr;
end
numBlocksTry_ = 10;
%-------------------------------
% The unadorned SDU0 (or X0Err) are from the svd. SDU0_ and X0Err_ are from
% subequent iterations of minimizing the collapsed objective.
SDU0Init = zeros(length(numBlocksTry_),MC); SDU0_ = zeros(length(numBlocksTry_),MC);
SDU0CllpsInit = zeros(length(numBlocksTry_),MC);  SDU0Cllps_ = zeros(length(numBlocksTry_),MC); 
SDU0Perm = zeros(length(numBlocksTry_),MC); SDU0Perm_ = zeros(length(numBlocksTry_),MC);
X0Err = zeros(length(numBlocksTry_),MC); X0Err_ = zeros(length(numBlocksTry_),MC);
X0CllpsErr = zeros(length(numBlocksTry_),MC); X0CllpsErr_ = zeros(length(numBlocksTry_),MC);
X0PermErr = zeros(length(numBlocksTry_),MC);  X0PermErr_ = zeros(length(numBlocksTry_),MC);
%fill = "sum"; %fill = "median"; %fill = "both";
fill = "one";
%fill = "mean";
%-----------------------
if real % ground truth matrix read only once if the data are real
    idxFlag = 1; % idxFlag = 1 means generate indices for synthetic data  
    [Xzeros, rowIdx, colIdx, Xcol, Xrow,idx] = processMatrix(X, n, q, p,real,idxFlag,'idxPlaceholder');
end
for i = 1 : length(numBlocksTry_)
    numBlocksTry = numBlocksTry_(i);
    for mc = 1 : MC
        if real == 0 % new matrix generated at every mc run only if synthetic data, for real data MC = 1 always
            idxFlag = 1; 
            [Xzeros, rowIdx, colIdx, Xcol, Xrow,idx] = processMatrix(X, n, q, p,real,idxFlag,'idxPlaceHolder');    
        end
        % Collapsed
        [XzerosPerm, XzerosCllps, r_All] = processBlocks(rowIdx, Xcol, Xzeros, q, numBlocksTry, same, fill); 
        % Method 1. U0 from unpermuted
        %[U0, S0, V0] = svds(Xzeros,r);
        %SDU0Init(i,mc) = norm(Ustr - U0*(U0'*Ustr));        
        %X0 = U0*S0*V0';
        %X0Err(i,mc) = norm(X0(idx) - Xzeros(idx))/norm(Xzeros(idx));
        % U0 from unpermuted plus iterations
        %P_updt = 0;
        %[SDVals,Errs] = altGDMin_LRMC(n,q,r,'r_PlaceHolder',p,U0, ...
        %                              Ustr,T, ...
        %                              rowIdx,Xcol,colIdx,Xrow, X0,idx,Xzeros,real, P_updt);        
        %SDU0_(i,mc) = SDVals(end);   
        %X0Err_(i,mc) = Errs(end);  
        %if real
        %    Ustr = U0;
        %end
        %{
        %------------------------------------
        % Method 2. U0 from cllps
        [U0Cllps, S0Cllps, V0Cllps] = svds(XzerosCllps,r);
        U0Cllps = U0Cllps(:,1:r); S0Cllps = S0Cllps(1:r,1:r); V0Cllps = V0Cllps(:,1:r);
        X0Cllps = U0Cllps*S0Cllps*V0Cllps';
        X0CllpsErr(i,mc) = norm(X0Cllps(idx) - Xzeros(idx))/norm(Xzeros(idx));
        SDU0Cllps(i,mc) = norm(Ustr - U0Cllps*(U0Cllps'*Ustr));
        % U0 from cllps matrix plus iterations
        idxFlag = 0; 
        [XzerosCllps, rowIdxCllps, colIdxCllps, XcolCllps, XrowCllps,~] = processMatrix(XzerosCllps, n, q, p,real,idxFlag,idx);
        [SDVals,Errs] = altMinwithP_LRMC(n,q, r, U0Cllps, Ustr, T, ...
                                       rowIdxCllps, XcolCllps, colIdxCllps, XrowCllps, ...
                                       X0Cllps,idx,Xzeros,real);
        X0CllpsErr_(i,mc) = Errs(end);
        SDU0Cllps_(i,mc) = SDVals(end);  
        % Method 3. U0 from permuted matrix
        [U0Perm, S0Perm, V0Perm] = svds(XzerosPerm,r);
        U0Perm = U0Perm(:,1:r);        S0Perm = S0Perm(1:r,1:r);        V0Perm = V0Perm(:,1:r);
        X0Perm = U0Perm*S0Perm*V0Perm';        
        X0PermErr(i,mc) = norm(X0Perm(idx) - Xzeros(idx))/norm(Xzeros(idx));
        SDU0Perm(i,mc) = norm(Ustr - U0Perm*(U0Perm'*Ustr));    
        % U0 from permuted matrix plus iterations        
        idxFlag = 0;
        [XzerosPerm, rowIdxPerm, colIdxPerm, XcolPerm, XrowPerm,~] = processMatrix(XzerosPerm, n, q, p,real,idxFlag,idx);        
        [SDVals,Errs] = altMinInit_LRMC(n,q, r, U0Perm, Ustr, T, ...
                                       rowIdxPerm, XcolPerm, colIdxPerm, XrowPerm, ...
                                       X0Perm,idx,Xzeros,real);
        X0PermErr_(i,mc) = Errs(end);    
        SDU0Perm_(i,mc) = SDVals(end);
        %}
        [U0Cllps, S0Cllps, V0Cllps] = svds(XzerosCllps,r);
        X0Cllps = U0Cllps*S0Cllps*V0Cllps';
        X0CllpsErr(i,mc) = norm(X0Cllps(idx) - Xzeros(idx))/norm(Xzeros(idx));
        [Uproj,~,~] = qr(U0Cllps,'econ');
        idxFlag = 0;
        [XzerosPerm, rowIdxPerm, colIdxPerm, XcolPerm, XrowPerm,~] = processMatrix(XzerosPerm, n, q, p,real,idxFlag,idx);
        T_in = 100;
        [SDVals,objVals,U0Cllps] = minCllpsInit(n,q,r,r_All,U0Cllps, ...
                                        Ustr,Bstr,T_init, ...
                                        rowIdxPerm,XcolPerm, 'zero','zero','zero',0,T_in);
        SDU0CllpsInit(i,mc) = norm(Ustr - U0Cllps*(U0Cllps'*Ustr));
        P_updt = 1;
        altMin = 1;
        [SDVals,Errs] = altGDMin_LRMC(n,q,r,r_All,p,U0Cllps, ...
                                      Ustr,Bstr,T_AltMin, ...
                                      rowIdxPerm,XcolPerm,colIdx,Xrow, X0Cllps,idx,XzerosPerm,real,P_updt,altMin);                
        SDU0Cllps_(i,mc) = SDVals(end);        
        mc
    end
end
SDU0Init = sum(SDU0Init,2)/MC; SDU0_ = sum(SDU0_,2)/MC;
SDU0CllpsInit = sum(SDU0CllpsInit,2)/MC; SDU0Cllps_ = sum(SDU0Cllps_,2)/MC;
SDU0Perm = sum(SDU0Perm,2)/MC; SDU0Perm_ = sum(SDU0Perm_,2)/MC;
%---------------------------------
X0Err = sum(X0Err,2)/MC; X0Err_ = sum(X0Err_,2)/MC;
X0CllpsErr = sum(X0CllpsErr,2)/MC; X0CllpsErr_ = sum(X0CllpsErr_,2)/MC;
X0PermErr = sum(X0PermErr,2)/MC; X0PermErr_ = sum(X0PermErr_,2)/MC;
%---
%SDU0_'
%SDU0Cllps_'
%SDU0Perm_'
%SDInit
X0Err_'
X0CllpsErr_'
X0PermErr_'
if real == 0
    plotRsltsLRMC(SDU0_, 0, SDU0Cllps_, SDU0Perm_,n,q,r,p,numBlocksTry_,MC,same,fill,real,T_init);
end
n,q
toc
function [Xzeros, rowIdx, colIdx, Xcol, Xrow,idx] = processMatrix(X, n, q, p,real,idxFlag,idx)
    if real
        idx = find(X > 0);
    elseif idxFlag      % Randomly select indices based on probability p
        idx = randperm(n * q);
        idx = idx(1:round(p * n * q));
    end
    % Convert linear indices to subscripts
    [row, col] = ind2sub([n, q], idx);
    % Instantiate Y = X_Omega (Observed entries)
    Xzeros = zeros(n, q);
    Xzeros(idx) = X(idx);
    % Initialize cell arrays
    rowIdx = cell(q, 1); 
    colIdx = cell(n, 1);
    Xcol = cell(q, 1); 
    Xrow = cell(n, 1);
    % Parallel processing
    for j = 1:q
        rowIdx{j} = row(col == j);
        Xcol{j} = X(rowIdx{j}, j);
    end
    for j = 1:n
        colIdx{j} = col(row == j);
        Xrow{j} = X(j, colIdx{j})';        
    end
end
%--------------------------
function [XzerosPerm, XzerosCllps,r_All] = processBlocks(rowIdx, Xcol, Xzeros, q, numBlocksTry, same, fill)
    % Initialize output variables
    XzerosPerm = Xzeros;
    XzerosCllps = Xzeros;
    r_All = cell(q,1);
    for k = 1 : q
        l_k = length(rowIdx{k});
        numBlocks = min(numBlocksTry, l_k);
        r_ = floor(l_k / numBlocks) * ones(numBlocks, 1);      
        rem = mod(l_k,numBlocks);
        if rem > 0
            r_(1:rem) = r_(1:rem) + 1; %r_ is a vector of length numBlocks, 
            %r_(end) = r_(end) + l_k - floor(l_k / numBlocks) * numBlocks;
        end
        r_All{k} = r_;
        if ~same
            pi_map = get_permutation_r(l_k, r_);
        else
            pi_map = 1:l_k; % Identity mapping if 'same' is true
        end        
        rowIdxPerm_k = rowIdx{k}(pi_map);
        XcolPerm_k = Xzeros(rowIdxPerm_k, k);
        XzerosPerm(rowIdxPerm_k, k) = Xcol{k};
        start  = 1;
        for s = 1 : numBlocks
            stop = start + r_(s) - 1;
            if fill == "mean"
                avg = sum(XcolPerm_k(start:stop)) / r_(s);
                XzerosCllps(rowIdxPerm_k(start:stop), k) = avg;                
            elseif fill == "median"
                XzerosCllps(rowIdxPerm_k(start:stop), k) = median(XcolPerm_k(start:stop));
            elseif fill == "one"
                XzerosCllps(rowIdxPerm_k(start), k) = sum(XcolPerm_k(start:stop));                
            elseif fill == "sum"
                XzerosCllps(rowIdxPerm_k(start:stop), k) = sum(XcolPerm_k(start:stop));                
            else
                XzerosCllps(rowIdxPerm_k(start:stop), k) = sum(XcolPerm_k(start:stop)) / r_(s);
                XzerosCllps(rowIdxPerm_k(floor(start + (stop - start) / 2):stop), k) = median(XcolPerm_k(start:stop));
            end
            start = start + r_(s);
        end
    end
end
