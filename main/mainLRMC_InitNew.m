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
r = 5;
n = 1000;
q = 1000;
p = 0.1;
same = 0;
numBlocks = n/10;   %effectively, m_new = numBlocks
r_ = ones(1,numBlocks)*(n/numBlocks);
MC = 1;
% generate rank-r X*
Ustr = orth(randn(n,r));
Bstr = randn(r,q);
X  = Ustr*Bstr;
% generate q matrices A_k of size m x n, m << n
%SDVals_UnPerm = zeros(MC,T+1);
%SDVals_sLcl = zeros(MC,T+1);
%SDVals_Perm = zeros(MC,T+1);
%------------------------
rowIdxPerm = cell(q,1);
XcolPerm = cell(q, 1); 
XzerosPerm = zeros(n,q);
XzerosCllps = zeros(n,q);
%XzerosCllps_ = zeros(numBlocks,q);
AkCllps_ = cell(q,1);
MCllpsSpctrl = zeros(n,q);
%-----------------------
for mc = 1 : MC
    [Xzeros, rowIdx, colIdx, Xcol, Xrow] = processMatrix(X, n, q, p);
    for k = 1 : q
        if ~same
            pi_map = get_permutation_r(n,r_);
        end
        for s = 1 : numBlocks
            start = sum(r_(1:s)) - r_(s) + 1;
            stop = start + r_(s) - 1;
            % --- Either 
            idx = start:stop;
            idx = idx(Xzeros(start:stop,k) ~= 0);
            XzerosCllps(idx,k) = sum(Xzeros(start:stop,k))/r_(s);
            % --- Or            
            %XzerosCllps_(s,k) = sum(XcolPerm{k}(start:stop));
            %AkCllps_{k}(s,start:stop) = ones(1,r_(s));
        end
    end
    [U0,~,~,] = svd(Xzeros,"econ");
    U0 = U0(:,1:r);
    SDU0 = norm(Ustr - U0*(U0'*Ustr))
    %------------------------------------------
    %[U0CllpsSpctrl,~,~,] = svd(MCllpsSpctrl,"econ");
    %U0CllpsSpctrl = U0CllpsSpctrl(:,1:r);
    %SDU0CllpsSpctrl = norm(Ustr - U0CllpsSpctrl*(U0CllpsSpctrl'*Ustr))
    %-----------------------------------
    %norm(Xzeros - XzerosCllps,'fro')
    [U0Cllps,~,~] = svd(XzerosCllps,"econ");
    U0Cllps = U0Cllps(:,1:r);
    SDU0Cllps = norm(Ustr - U0Cllps*(U0Cllps'*Ustr))
    %---------------------------------------
    %[U0Perm,~,~] = svd(XzerosPerm,"econ");
    %U0Perm = U0Perm(:,1:r);
    %SDU0Perm = norm(Ustr - U0Perm*(U0Perm'*Ustr))    
    mc
end
function [Xzeros, rowIdx, colIdx, Xcol, Xrow] = processMatrix(X, n, q, p)
    % Randomly select indices based on probability p
    idx = randperm(n * q);
    idx = idx(1:round(p * n * q));
    idxC = setdiff(1:n * q, idx);

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