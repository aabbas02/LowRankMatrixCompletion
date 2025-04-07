function [SDVals,objVals] = altMinwithP_LRMC(n,q,r,Uinit, ...
                                            Ustr,T, ...
                                            rowIdx,Xcol,colIdx,Xrow, Xhat0,idx,Xzeros,real)
    % IMPORTANT. DOES NOT WORK WITH P-update
    % NOTE: P update does not work. Set Pupdt = 0 manually below.
    % NOTE: P update does not work. Set Pupdt = 0 manually below.
    % NOTE: P update does not work. Set Pupdt = 0 manually below.
    % Explain why P update does not work.
    % altGDMinwithP_LRMC works.
    Pupdt = 0;
    U = Uinit;
    if real == 0
        SDVals = zeros(T+1,1);
        SDVals(1) = norm(Ustr - U*(U'*Ustr));
    end
    objVals = zeros(T+1,1);
    objVals(1) = norm(Xzeros(idx) - Xhat0(idx))/norm(Xzeros(idx));   
    V = zeros(r,q);    
    for i = 1 : T
        for j = 1 : q
            if Pupdt == 1 && i > 1
                % P - update 
                xHat = U(rowIdx{j},:)*V(:,j);
                [~,idx1] = sort(xHat);
                [~,idx2] = sort(Xcol{j});
                rowIdx{j}(idx2) = rowIdx{j}(idx1);
                % V - update
                V(:,j) = U(rowIdx{j},:)\Xcol{j};
            else
                V(:,j) = U(rowIdx{j},:)\Xcol{j};
            end
        end
        % U update 
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
