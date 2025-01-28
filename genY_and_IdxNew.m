function [Y,Ycol,rowIdx,Ycol_,rowsJ] = genY_and_IdxNew(Xstar, p, numWrkrs)
    [n, q] = size(Xstar);
    Y = zeros(n, q);  % Initialize Y as a zero matrix
    
	idx = randperm(n*q);
	idx = idx(1:round(p*n*q));
	%idxC = setdiff(1:n*q,idx);
	Y(idx) = Xstar(idx);
	[row,col] = ind2sub([n,q],idx);
	% Make cell arrays for row and column indices
	rowIdx = cell(q,1); %colIdx = cell(n,1);
	Ycol = cell(q,1); %Xrow = cell(n,1);
    for k = 1 : q
		rowIdx{k} = row(col==k);
		Ycol{k} =  Xstar(rowIdx{k},k);
		%if k <= n
		%	colIdx{k} = col(row==k);
		%	Xrow{k} = X(k,colIdx{k})';
		%end
		Y(rowIdx{k},k) = Xstar(rowIdx{k},k);
    end
    Ycol_ = cell(numWrkrs,q/numWrkrs);
    rowsJ = reshape(rowIdx,[q/numWrkrs,numWrkrs]);
    rowsJ = rowsJ';
    for j = 1 : numWrkrs
        offset = (j-1)*q/numWrkrs;
        for k = 1 : q/numWrkrs
            Ycol_{j,k} = Ycol{offset+k};
        end
    end
end
