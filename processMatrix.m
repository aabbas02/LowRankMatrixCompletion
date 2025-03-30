function [Xzeros, rowIdx, colIdx, Xcol, Xrow,idx, Xcol_,rowsJ] = processMatrix(X, p,real,idxFlag,L)
    n = size(X,1); q = size(X,2);
    if real
        idx = find(X > 0);
    elseif idxFlag      % Randomly select indices based on probability p
        idx = randperm(n * q);
        idx = idx(1:round(p * n * q));
    end
    % Convert linear indices to subscripts
    [row, col] = ind2sub([n, q], idx);
 
    % Instantiate Y = X_Omega (Observed entries)clclc
    Xzeros = zeros(n, q);
    Xzeros(idx) = X(idx);
    q = L*floor(q/L);
    Xzeros = Xzeros(:, 1:q);
    idx = find(Xzeros > 0);
    % Initialize cell arrays
    rowIdx = cell(q, 1); 
    colIdx = cell(n, 1);
    Xcol = cell(q, 1); 
    Xrow = cell(n, 1);
    % Parallel processing
    for j = 1 : q
        rowIdx{j} = row(col == j);
        Xcol{j} = X(rowIdx{j}, j);
        if mod(j,1000) == 0
            j
        end
    end
    disp("cols done")
    for j = 1 : n
        colIdx{j} = col(row == j);
        Xrow{j} = X(j, colIdx{j})';       
        if mod(j,1000) == 0
            j
        end        
    end
    disp("rows done")
    Xcol_ = cell(L,q/L);  
    rowsJ = reshape(rowIdx,[q/L,L]);
    rowsJ = rowsJ';    
    for j = 1 : L
        offset = (j-1)*q/L;
        for k = 1 : q/L
            Xcol_{j,k} = Xcol{offset+k};
        end
    end    
end