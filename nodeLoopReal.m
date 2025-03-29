function [grad_] = nodeLoopReal(U,q,n,r,rowIdx,Xcol,numWrkrs,rowsJ,Xcol_)
	% rowIdx is a cell array of size (1 x q)
	% colIdx is a cell array of size (1 x n)
	% Xcol is a cell array of size (1 x q)
	%%% rows are row indices of the observed entries, an array (NOT CELL ARRAY) of size (1 x number of observed entries)
	%%% cols are the corresponding column indices of the observed entries, an array (NOT CELL ARRAY) of size (1 x number of observed entries)
	%Grad_U = zeros(n,r);
    grad_ = zeros(n,r,numWrkrs);  
    for j = 1 : numWrkrs
        tmp = zeros(r,q/numWrkrs);
        diff = [];
        rows = [];
        cols = [];
        rowIdxj = rowsJ(j,:); 
                            % rowIdxj should be cell aray of size  1 x q/4.
                            % alternatively, rows_{j} can be a cell
                            % array of size 4 x 1,
        XcolJ = Xcol_(j,:);
        strt = 0;
        for k = 1 : q/numWrkrs
           rowIdx_jk = rowIdxj{k};
           numRows = length(rowIdx_jk);
           tmp(:,k) = U(rowIdx_jk,:)\XcolJ{k};
           diff(strt + 1 : strt + numRows) = U(rowIdx_jk,:)*tmp(:,k) - XcolJ{k}; 
           rows(strt + 1 : strt + numRows) = rowIdx_jk;
           cols(strt + 1 : strt + numRows) = k;
           strt = strt + numRows;
        end
        diff = diff(1:strt);
        rows = rows(1:strt);
        cols = cols(1:strt);
        M = sparse(rows,cols,diff,n,q/numWrkrs);
        grad_(:,:,j) = M*tmp'; 
    end
