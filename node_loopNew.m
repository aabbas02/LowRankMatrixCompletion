function Grad_U = node_loopNew(U,q,n,r,rowIdx,Xcol)
	% rowIdx is a cell array of size (1 x q)
	% colIdx is a cell array of size (1 x n)
	% Xcol is a cell array of size (1 x q)
	%%% rows are row indices of the observed entries, an array (NOT CELL ARRAY) of size (1 x number of observed entries)
	%%% cols are the corresponding column indices of the observed entries, an array (NOT CELL ARRAY) of size (1 x number of observed entries)
	%Grad_U = zeros(n,r);
	tmp = zeros(r,q);
	diff = [];
    rows = [];
    cols = [];
	strt = 0;
	for k = 1 : q
	   rowIdx_jk = rowIdx{k};
	   numRows = length(rowIdx_jk);
	   tmp(:,k) = U(rowIdx_jk,:)\Xcol{k};
	   diff(strt + 1 : strt + numRows) = U(rowIdx_jk,:)*tmp(:,k) - Xcol{k};
       rows(strt + 1 : strt + numRows) = rowIdx_jk;
       cols(strt + 1 : strt + numRows) = k;
	   strt = strt + numRows;
	end
	diff = diff(1:strt);
    rows = rows(1:strt);
    cols = cols(1:strt);
	M = sparse(rows,cols,diff,n,q);
	Grad_U = M*tmp'; 