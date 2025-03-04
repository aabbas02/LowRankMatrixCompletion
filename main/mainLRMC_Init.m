clc
close all
clear all
dir = pwd;
% For linux, replace '\' with '/'
cd ..
%addpath(genpath('.\functionsMtrxSnsng'));
addpath(genpath('.\utils'));
cd(dir)    
tic
%---------------------------------
r = 5;
real = 0;
T = 0;
n = 1000;
q = 1000;
p = 0.2;
same = 0;
MC = 5;
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
numBlocksTry_ = 10:10:100;
%-------------------------------
SDU0 = zeros(length(numBlocksTry_),MC);
SDU0_ = zeros(length(numBlocksTry_),MC);
SDU0Cllps = zeros(length(numBlocksTry_),MC); 
SDU0Cllps_ = zeros(length(numBlocksTry_),MC); 
SDU0Perm = zeros(length(numBlocksTry_),MC);
SDU0Perm_ = zeros(length(numBlocksTry_),MC);
X0Err = zeros(length(numBlocksTry_),MC);
X0Err_ = zeros(length(numBlocksTry_),MC);
X0CllpsErr = zeros(length(numBlocksTry_),MC);
X0CllpsErr_ = zeros(length(numBlocksTry_),MC);
X0PermErr = zeros(length(numBlocksTry_),MC); 
X0PermErr_ = zeros(length(numBlocksTry_),MC);
fill = "mean";
%fill = "median";
%fill = "both";
%-----------------------
if real % ground truth matrix read only once if the data are real
    idxFlag = 1;
    [Xzeros, rowIdx, colIdx, Xcol, Xrow,idx] = processMatrix(X, n, q, p,real,idxFlag,0);
end
for i = 1 : length(numBlocksTry_)
    numBlocksTry = numBlocksTry_(i);
    for mc = 1 : MC
        if real == 0 % new matrix generated at every mc run only if synthetic data
            idxFlag = 1;
            [Xzeros, rowIdx, colIdx, Xcol, Xrow,idx] = processMatrix(X, n, q, p,real,idxFlag,0);    
        end
        % U0 from unpermuted
        [U0, S0, V0] = svds(Xzeros,r);
        U0 = U0(:,1:r);
        S0 = S0(1:r,1:r);
        V0 = V0(:,1:r);
        X0 = U0*S0*V0';
        X0Err(i,mc) = norm(X0(idx) - Xzeros(idx))/norm(Xzeros(idx));
        % U0 from unpermuted plus iterations
        Pupdt = 0;
        [SDVals,Errs] = altMinInit(n,q, r, U0, Ustr, T, rowIdx, Xcol, colIdx, Xrow, X0,idx,Xzeros,real,Pupdt);
        SDU0_(i,mc) = SDVals(end);        
        X0Err_(i,mc) = Errs(end);  
        if real
            Ustr = U0;
        end
        SDU0(i,mc) = norm(Ustr - U0*(U0'*Ustr));
        % get permuted and collapsed matrices -- used below
        [XzerosPerm, XzerosCllps] = processBlocks(rowIdx, Xcol, Xzeros, q, numBlocksTry, same, fill);    
        %------------------------------------
        % ----------- U0 from cllps
        [U0Cllps, S0Cllps, V0Cllps] = svds(XzerosCllps,r);
        U0Cllps = U0Cllps(:,1:r);
        S0Cllps = S0Cllps(1:r,1:r);
        V0Cllps = V0Cllps(:,1:r);
        X0Cllps = U0Cllps*S0Cllps*V0Cllps';
        X0CllpsErr(i,mc) = norm(X0Cllps(idx) - Xzeros(idx))/norm(Xzeros(idx));
        %SDU0Cllps(i,mc) = norm(Ustr - U0Cllps*(U0Cllps'*Ustr));
        % U0 from cllps matrix plus iterations
        idxFlag = 0;
        Pupdt = 0;
        [XzerosCllps, rowIdxCllps, colIdxCllps, XcolCllps, XrowCllps,~] = processMatrix(XzerosCllps, n, q, p,real,idxFlag,idx);
        [SDVals,Errs] = altMinInit(n,q, r, U0Cllps, Ustr, T, ...
                                       rowIdxCllps, XcolCllps, colIdxCllps, XrowCllps, ...
                                       X0Cllps,idx,Xzeros,real,Pupdt);
        SDU0Cllps_(i,mc) = SDVals(end);
        X0CllpsErr_(i,mc) = Errs(end);
        % ------------ U0 from permuted matrix
        [U0Perm, S0Perm, V0Perm] = svds(XzerosPerm,r);
        U0Perm = U0Perm(:,1:r);
        S0Perm = S0Perm(1:r,1:r);
        V0Perm = V0Perm(:,1:r);
        X0Perm = U0Perm*S0Perm*V0Perm';        
        X0PermErr(i,mc) = norm(X0Perm(idx) - Xzeros(idx))/norm(Xzeros(idx));
        %SDU0Perm(i,mc) = norm(Ustr - U0Perm*(U0Perm'*Ustr));    
        % U0 from permuted matrix plus iterations        
        idxFlag = 0;
        Pupdt = 0;
        [XzerosPerm, rowIdxPerm, colIdxPerm, XcolPerm, XrowPerm,~] = processMatrix(XzerosPerm, n, q, p,real,idxFlag,idx);        
        [SDVals,Errs] = altMinInit(n,q, r, U0Perm, Ustr, T, ...
                                       rowIdxPerm, XcolPerm, colIdxPerm, XrowPerm, ...
                                       X0Perm,idx,Xzeros,real,Pupdt);
        SDU0Perm_(i,mc) = SDVals(end);
        X0PermErr_(i,mc) = Errs(end);
        
        mc
    end
end
SDU0 = sum(SDU0,2)/MC;
SDU0_ = sum(SDU0_,2)/MC;
SDU0Cllps = sum(SDU0Cllps,2)/MC;
SDU0Cllps_ = sum(SDU0Cllps_,2)/MC;
SDU0Perm = sum(SDU0Perm,2)/MC;
SDU0Perm_ = sum(SDU0Perm_,2)/MC;
%---
SDU0_'
SDU0Cllps_'
SDU0Perm_'

X0Err_'
X0CllpsErr_'
X0PermErr_'
if real == 0
    plotRslts(SDU0_, 0, SDU0Cllps_, SDU0Perm_,n,q,r,p,numBlocksTry_,MC,same,fill,real,T);
end
n,q
toc
function [XzerosPerm, XzerosCllps] = processBlocks(rowIdx, Xcol, Xzeros, q, numBlocksTry, same, fill)
    % Initialize output variables
    XzerosPerm = Xzeros;
    rowIdxPerm = cell(q,1);
    XcolPerm = cell(q, 1); 
    XzerosCllps = Xzeros;

    for k = 1 : q
        l_k = length(rowIdx{k});
        numBlocks = min(numBlocksTry, l_k);
        r_ = floor(l_k / numBlocks) * ones(numBlocks, 1);        
        if mod(l_k, numBlocks) > 0
            r_(end) = r_(end) + l_k - floor(l_k / numBlocks) * numBlocks;
        end        
        if ~same
            pi_map = get_permutation_r(l_k, r_);
        else
            pi_map = 1:l_k; % Identity mapping if 'same' is true
        end        
        rowIdxPerm{k} = rowIdx{k}(pi_map);
        XcolPerm{k} = Xzeros(rowIdxPerm{k}, k);
        XzerosPerm(rowIdxPerm{k}, k) = Xcol{k};
        for s = 1:numBlocks
            start = sum(r_(1:s)) - r_(s) + 1;
            stop = start + r_(s) - 1;
            if fill == "mean"
                avg = sum(XcolPerm{k}(start:stop)) / r_(s);
                XzerosCllps(rowIdxPerm{k}(start:stop), k) = avg;
            elseif fill == "median"
                XzerosCllps(rowIdxPerm{k}(start:stop), k) = median(XcolPerm{k}(start:stop));
            else
                XzerosCllps(rowIdxPerm{k}(start:stop), k) = sum(XcolPerm{k}(start:stop)) / r_(s);
                XzerosCllps(rowIdxPerm{k}(floor(start + (stop - start) / 2):stop), k) = median(XcolPerm{k}(start:stop));
            end
        end
    end

end


function plotRslts(SDU0, SDU0Cllps, SDU0Cllps_, SDU0Perm,n,q,r,p,numBlocks_,MC,same,fill,real,T)
    figure;
    hold on
    if real == 0
    plot(numBlocks_,SDU0, ...
        'DisplayName', 'SDVals Unperm', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    end
    %plot(numBlocks_,SDU0Cllps, ...
    %    'DisplayName', 'SDVals Cllps', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    plot(numBlocks_,SDU0Cllps_, ...
        'DisplayName', 'SDVals Cllps Mean Only', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    plot(numBlocks_,SDU0Perm, ...
        'DisplayName', 'SDVals Naive', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    grid("on")
    xticks(numBlocks_);
    %-------------------------------
    title("LRMC. n = " + n + ", q = " + q +...
          ", r = " + r + ", p = " + p  +  ", MC = " + MC + ", same = " + same + ", fill = " + fill + ", T = " + T,...
           'Interpreter', 'Latex', 'FontSize',11)
    %--------------------------------
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    if real == 0
        ylabel("$SD(U^{(T)},U^*)$","FontSize",14,'Interpreter','Latex')
    else
        ylabel("Initialization Error", "FontSize",11,"Interpreter","Latex")
    end
    xlabel('number of blocks', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['LRMC_Init_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_p_', num2str(p),'_numBlocks_', num2str(max(numBlocks_)), '_same_',num2str(same),'_fill_',num2str(fill), '_T_',num2str(T)];
    
    savefig([stringTitle, '.fig']);
end
%---
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
%---
function [Ustr,X,p] = getMovieLens(r)
    A = readmatrix("ratings100K.xlsx");
    %--------------------
    %load("ratings1M.mat");
    %--------------------
    %load("ratings10M.mat")
    n = max(A(:,1));
    q = max(A(:,2));
    X = zeros(n,q);
    num = 0;
    disp(num)
    for k = 1 : size(A,1)
        i = A(k,1);
        j = A(k,2);
        X(i,j) = A(k,3);
        num = num + 1;
        %if mod(k,100000) == 0
        %    size(A,1),k
        %end
    end
    p = num/(n*q);
    %X = X';
    %X = X(1:10000,1:10000);
    [Ustr,~,~] = svds(X,r);
    Ustr = Ustr(:,1:r);
end
