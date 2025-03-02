function [SDVals,objVals] = altMinWithP_LRMC(n,q,r,Uinit, ...
                                        Ustr,T, ...
                                        rowIdx,Xcol,colIdx,Xrow, Xhat0,idx,Xzeros,real,B_init)
    % Given U^{(0)}, B^{(0)}
    % Update P_j for all j in [q]. This is to to update rowIdx{j}
    % Update b_j for all j in [q]. This is a least-squares problem with
    % P_jU as the sensing matrix.
    % Update U - this is the tricky part because
    % The first row of U does not depend on the observed entries of X* in
    % row 1 anymore ( that is after the P update).

    V = B_init;
    U = Uinit;
    if real == 0
        SDVals = zeros(T+1,1);
        SDVals(1) = norm(Ustr - U*(U'*Ustr));
    end

    objVals = zeros(T+1,1);
    objVals(1) = norm(Xzeros(idx) - Xhat0(idx))/norm(Xzeros(idx));   
    for i = 1 : T
        % P update
        for j = 1 : q
            xHat = U(rowIdx{j},:)*V(:,j);
            [~,idx1] = sort(xHat);
            [~,idx2] = sort(Xcol{j}); % makre sure this does not do an in-place sort
            %V(:,j) = U(rowIdx{j}(idx1),:)\Xcol{j}(idx2);                
            rowIdx{j}(idx2) = rowIdx{j}(idx1);
        end
        % V update
        for j = 1 : q
            V(:,j) = U(rowIdx{j},:)\Xcol{j};
        end
        % U update 
        for j = 1 : n
			U(j,:) = V(:, colIdx{j})'\???; % This right hand side should change 
        end
        if real == 0
            [Uproj,~,~] = qr(U,'econ');
            SDVals(i + 1) = norm(Ustr - Uproj*(Uproj'*Ustr));  
        end
        Xhat = U*V;
        objVals(i+1) = norm(Xzeros(idx) - Xhat(idx))/norm(Xzeros(idx));
    end
