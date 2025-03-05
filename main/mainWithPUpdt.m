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
T_init = 100;
T = 500;
n = 1000;
q = 1000;
p = 0.05;
same = 0;
MC = 25;
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
numBlocksTry_ = [20];
SDValsAltGDMin = zeros(MC,T+1);
%-------------------------------
fill = "mean";
%fill = "median";
%fill = "both";
%-----------------------
if real % ground truth matrix read only once if the data ar real
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
        [XzerosPerm, XzerosCllps,r_] = processBlocks(rowIdx, Xcol, Xzeros, q, numBlocksTry, same, fill);  
        idxFlag = 0;
        [XzerosCllps, rowIdxCllps, colIdxCllps, XcolCllps, XrowCllps,~] = processMatrix(XzerosCllps, n, q, p,real,idxFlag,idx);
        [U0Cllps, S0Cllps, V0Cllps] = svds(XzerosCllps,r);        
        Pupdt = 0;
        [SDVals,~,U_init,B_init] = altMinInit(n,q, r, U0Cllps, Ustr, T_init, ...
                                       rowIdxCllps, XcolCllps, colIdxCllps, XrowCllps, ...
                                       U0Cllps*S0Cllps*V0Cllps',idx,XzerosCllps,real,Pupdt);
        SDVals
        
        X_init = U_init*B_init;
        [~, rowIdxPerm, colIdxPerm, XcolPerm, XrowPerm,~] = processMatrix(XzerosPerm, n, q, p,real,idxFlag,idx);
        [SDValsAltGDMin(mc,:),Errs] = altGDMinWithP_LRMC(n,q, r,r_,p, U_init, Ustr, T, ...
                                    rowIdxPerm, XcolPerm, colIdxPerm, XrowPerm, ...
                                    X_init,idx,XzerosPerm,real,B_init);        
        mc
    end
end
SDValsAltGDMin = sum(SDValsAltGDMin,1)/MC;
if real == 0
    plotRslts(SDValsAltGDMin,n,q,r,p,numBlocksTry_(1),MC,same,fill,real,T,T_init);
end

function [XzerosPerm, XzerosCllps,r_cell] = processBlocks(rowIdx, Xcol, Xzeros, q, numBlocksTry, same, fill)
    % Initialize output variables
    XzerosPerm = Xzeros;
    rowIdxPerm = cell(q,1);
    XcolPerm = cell(q, 1); 
    r_cell = cell(q,1);
    XzerosCllps = Xzeros;
    for k = 1 : q

        l_k = length(rowIdx{k});
        numBlocks = min(numBlocksTry, l_k);
        r_ = floor(l_k / numBlocks) * ones(numBlocks, 1);
        
        if mod(l_k, numBlocks) > 0
            r_(end) = r_(end) + l_k - floor(l_k / numBlocks) * numBlocks;
        end
        r_cell{k} = r_;
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


function plotRslts(SDAltGDMin,n,q,r,p,numBlocks,MC,same,fill,real,T,T_init)
    figure;
    %hold on
    %if real == 0
    %plot(numBlocks_,SDU0, ...
    %    'DisplayName', 'SDVals Unperm', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    %end
    %plot(numBlocks_,SDU0Cllps, ...
    %    'DisplayName', 'SDVals Cllps', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    %plot(numBlocks_,SDU0Cllps_, ...
    %    'DisplayName', 'SDVals Cllps Mean Only', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    semilogy(SDAltGDMin, ...
        'DisplayName', 'SDVals AltMin', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    grid("on")
    %-------------------------------
    title("LRMC. n = " + n + ", q = " + q +...
          ", r = " + r + ", p = " + p  +  ", MC = " + MC + ", same = " + same + ", fill = " + fill + ", T = " + T + ", num Blocks = " + numBlocks,...
          'Interpreter', 'Latex', 'FontSize',11);
    %--------------------------------
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    if real == 0
        ylabel("$SD(U^{(t)},U^*)$","FontSize",14,'Interpreter','Latex')
    else
        ylabel("Initialization Error", "FontSize",11,"Interpreter","Latex")
    end
    xlabel('Iterations t', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['AltMinbyGD', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_p_', num2str(p),'_numBlocks_', num2str((numBlocks)), '_same_',num2str(same),'_fill_',num2str(fill), '_T_',num2str(T), '_T_init_',num2str(T_init)];
    
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
