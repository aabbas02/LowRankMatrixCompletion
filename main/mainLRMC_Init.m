clc
close all
clear all
dir = pwd;
% For linux, replace '\' with '/'
cd ..
%addpath(genpath('.\functionsMtrxSnsng'));
addpath(genpath('.\utils'));
cd(dir)    
%---------------------------------
r = 10;
n = 1000;
q = 1000;
p = 0.1;
same = 0;
r_ = ones(1,numBlocks)*(n/numBlocks);
MC = 1;
real = 1;
%------------------------
[Ustr,X] = getMovieLens(5);
n = size(X,1);
q = size(X,2);
% generate rank-r X*
%Ustr = orth(randn(n,r));
%Bstr = randn(r,q);
%X  = Ustr*Bstr;
%----------------------------------------
%numBlocksTry_ = [10,20,30,40,50,60,70,80,90,100];
numBlocksTry_ = [100,120,140,160,180,200]/2;
% generate rank-r X*
%Ustr = orth(randi(1000,[n,r]));
%Bstr = randn(r,q);
%X  = Ustr*Bstr;
%X = rand(n,q);
%[Ustr,Sstr,Vstr] = svd(X,'econ');
%Ustr = Ustr(:,1:r);
%Sstr = Sstr(1:r,1:r);
%Vstr = Vstr(:,1:r);
%X = Ustr*Sstr*Vstr';
SDU0 = zeros(length(numBlocksTry_),MC);
SDU0Cllps = zeros(length(numBlocksTry_),MC);
SDU0Perm = zeros(length(numBlocksTry_),MC);
fill = "mean";
fill = "both";
%------------------------
%-----------------------
for i = 1 : length(numBlocksTry_)
    numBlocksTry = numBlocksTry_(i);
    for mc = 1 : MC
        rowIdxPerm = cell(q,1);
        XcolPerm = cell(q, 1); 
        XzerosPerm = zeros(n,q);
        XzerosCllps = zeros(n,q);
        %if same
        %    pi_map = get_permutation_r(n,r_);
        %end
        [Xzeros, rowIdx, colIdx, Xcol, Xrow] = processMatrix(X, n, q, p,real);
        for k = 1 : q
            l_k = length(rowIdx{k});
            numBlocks = min(numBlocksTry,l_k);
            r_ = floor(l_k/numBlocks)*ones(numBlocks,1);
            if mod(l_k,numBlocks) > 0
                r_(end) = r_(end) + l_k - floor(l_k/numBlocks)*numBlocks;
            end
            if ~same
                pi_map = get_permutation_r(l_k,r_);
            end
            rowIdxPerm{k} = rowIdx{k}(pi_map);
            XcolPerm{k} = Xzeros(rowIdxPerm{k},k);
            XzerosPerm(rowIdxPerm{k},k) = Xcol{k};
            %AkCllps_{k} = zeros(numBlocks,n);
            for s = 1 : numBlocks
                start = sum(r_(1:s)) - r_(s) + 1;
                stop = start + r_(s) - 1;
                if fill == "mean"
                    % Either replace all entries by average
                    XzerosCllps(rowIdxPerm{k}(start:stop),k) = sum(XcolPerm{k}(start:stop))/r_(s) + 0*1e0*randn(length(start:stop),1);
                else
                    % Or replace half by average, half by median
                    XzerosCllps(rowIdxPerm{k}(start:stop),k) = sum(XcolPerm{k}(start:stop))/r_(s);
                    XzerosCllps(rowIdxPerm{k}(floor(start +(stop-start)/2):stop),k) = median(XcolPerm{k}(start:stop));
                end
                %AkCllps_{k}(s,start:stop) = ones(1,r_(s));
            end
        end
        [U0,~,~,] = svd(Xzeros,"econ");
        U0 = U0(:,1:r);
        SDU0(i,mc) = norm(Ustr - U0*(U0'*Ustr));
        %-----------------------------------
        %norm(Xzeros - XzerosCllps,'fro')
        [U0Cllps,~,~] = svd(XzerosCllps,"econ");
        U0Cllps = U0Cllps(:,1:r);
        SDU0Cllps(i,mc) = norm(Ustr - U0Cllps*(U0Cllps'*Ustr));
        %---------------------------------------
        [U0Perm,~,~] = svd(XzerosPerm,"econ");
        U0Perm = U0Perm(:,1:r);
        SDU0Perm(i,mc) = norm(Ustr - U0Perm*(U0Perm'*Ustr));    
        mc
    end
end
SDU0 = sum(SDU0,2)/MC;
SDU0Cllps = sum(SDU0Cllps,2)/MC;
SDU0Perm = sum(SDU0Perm,2)/MC;
%disp(C)
plotRslts(SDU0, SDU0Cllps, SDU0Perm,n,q,r,p,numBlocksTry_,MC,same,fill);

function plotRslts(SDU0, SDU0Cllps, SDU0Perm,n,q,r,p,numBlocks_,MC,same,fill)
    figure;
    plot(numBlocks_,SDU0, ...
        'DisplayName', 'SDVals U^{(0)}', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    hold on;
    plot(numBlocks_,SDU0Cllps, ...
        'DisplayName', 'SDVals U^{(0)}-Cllps', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    plot(numBlocks_,SDU0Perm, ...
        'DisplayName', 'SDVals U^{(0)}-Perm', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    grid("on")
    xticks(numBlocks_);
    %-------------------------------
    title("LRMC Init. n = " + n + ", q = " + q +...
          ", r = " + r + ", p = " + p  +  ", MC = " + MC + ", same = " + same + ", fill = " + fill,...
           'Interpreter', 'Latex', 'FontSize',11)
    %--------------------------------
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel("$SD(U^{(0)},U^*)$","FontSize",14,'Interpreter','Latex')
    xlabel('number of blocks', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['LRMC_Init_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_p_', num2str(p), '_same_',num2str(same),'_fill_',num2str(fill)];
    
    savefig([stringTitle, '.fig']);
end
function [Xzeros, rowIdx, colIdx, Xcol, Xrow] = processMatrix(X, n, q, p,real)
    % Randomly select indices based on probability p
    if real
        idx = find(X > 0);
    else
        idx = randperm(n * q);
        idx = idx(1:round(p * n * q));
        idxC = setdiff(1:n * q, idx);
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
    parfor j = 1:q
        rowIdx{j} = row(col == j);
        Xcol{j} = X(rowIdx{j}, j);
        if j <= n
            colIdx{j} = col(row == j);
            Xrow{j} = X(j, colIdx{j})';
        end
    end
end
%---
%plotRslts(SDVals_sLcl, SDVals_Perm, SDVals_UnPerm,n,q,r,m,numBlocks,MC,same);
function plotRsltsLRMCInit(numBlocks_, SDVals_,n,q,r,m,MC,same)
    figure;
    %SDVals_
    %SDVals_UnPerm = sum(SDVals_UnPerm,1)/MC;
    %SDVals_sLcl = sum(SDVals_sLcl,1)/MC;
    semilogy(numBlocks_,SDVals_, ...
        'DisplayName', 'SDVals-Spectral Init LRMC', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    hold on;
    %semilogy(SDVals_Perm, ...
    %    'DisplayName', 'SDVals (Naive)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 7);
    %
    %semilogy(SDVals_UnPerm, ...
    %    'DisplayName', 'SDVals (Unpermuted)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
    %grid on
    
    
    title("n = " + n + ", q = " + q +...
          ", r = " + r + ", m = " + m + ", num Blocks = " + numBlocks +  ", MC = " + MC + ", same = " + same, ...
           'Interpreter', 'Latex', 'FontSize',14)
    
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel("$SD(U,U^*)$","FontSize",14,'Interpreter','Latex')
    xlabel('Iterations (t)', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['MtrxSnsngPerm_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_m_', num2str(m), '_numBlocks_',num2str(numBlocks), '_same_',num2str(same)];
    
    savefig([stringTitle, '.fig']);
end
function [Ustr,X] = getMovieLens(r)
    A = readmatrix("ratings.xlsx");
    n = max(A(:,1));
    q = max(A(:,2));
    X = zeros(n,q);
    for k = 1 : size(A,1)
        i = A(k,1);
        j = A(k,2);
        X(i,j) = A(k,3);
    end
    [Ustr,~,~] = svd(X,"econ");
    Ustr = Ustr(:,1:r);
end
