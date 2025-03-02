function [SDVals,objVals] = altGDMin_LRMCPerm(n,q,r, ...
                                            Uinit, Ustr, T, ...
                                            rowIdx,Xcol,colIdx,Xrow, Xhat0, ...
                                            rowIdxCllps, XcolCllps, colIdxCllps, XrowCllps,X0Cllps,...
                                            idx,Xzeros,real, Pupdt)
    % 
    U = Uinit;
    if real == 0
        SDVals = zeros(T+1,1);
        SDVals(1) = norm(Ustr - U*(U'*Ustr));
    end
    objVals = zeros(T+1,1);
    objVals(1) = norm(Xzeros(idx) - Xhat0(idx))/norm(Xzeros(idx));   
    V = zeros(r,q);    
    for i = 1 : T
        % V update
        for j = 1 : q
            if i > 1 % update by full
                V(:,j) = U(rowIdx{j},:)\Xcol{j};
            else % update by collapsed
                V(:,j) = U(rowIdxCllps{j},:)\XcolCllps{j};
            end
        end
        % P update
        if Pupdt
            for k = 1 : q
                xHat = U(rowIdx{j},:)*V(:,j);
                [~,idx1] = sort(xHat);
                [~,idx2] = sort(Xcol{j});
                V(:,j) = U(rowIdx{j}(idx1),:)\Xcol{j}(idx2);                
                rowIdx{j}(idx2) = rowIdx{j}(idx1);                
            end
        end
        % U update by
        for j = 1 : n
			U(j,:) = V(:, colIdx{j})'\Xrow{j};
        end
        if real == 0
            [Uproj,~,~] = qr(U,'econ');
            SDVals(i + 1) = norm(Ustr - Uproj*(Uproj'*Ustr));  
        end
        Xhat = U*V;
        objVals(i+1) = norm(Xzeros(idx) - Xhat(idx))/norm(Xzeros(idx));
    end
