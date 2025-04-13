function [SDVals,objVals] = altGDMin_LRMC(n,q,r,r_All,p,Uinit, ...
                                        Ustr,Bstr,T, ...
                                        rowIdx,Xcol,colIdx,Xrow, Xhat0,idx,Xzeros,real,P_updt,altMin)
    % Given U^{(0)}, B^{(0)}
    % Update P_j for all j in [q]. This is to to update rowIdx{j}
    % Update b_j for all j in [q]. This is a least-squares problem with
    % P_jU as the sensing matrix.
    % Update U by gradient descent - Why is this easier than updating U by
    % least-squares?
    U = Uinit;
    if altMin
        T_in = 100;
    end
    if P_updt
        V_init = getVInit(r, q, r_All, rowIdx, U, Xcol);
    end
    gradU = zeros(n,r);
    [~,S] = svds(Xzeros/p,r);
    % Initialization of Uzero and Vzero after Projection 
    eta = 0.01/(S(1,1)^2*p); 
    if real == 0
        SDVals = zeros(T+1,1);
        SDVals(1) = norm(Ustr - U*(U'*Ustr));
    end
    objVals = zeros(T+1,1);
    objVals(1) = norm(Xzeros(idx) - Xhat0(idx))/norm(Xzeros(idx));   
    for i = 1 : T
        err = 0;
        % V update          
        if i == 1 && P_updt == 1 %collapsed least squares
            V = V_init; 
        else % full least squares
            for j = 1 : q
                V(:,j) = pinv(U(rowIdx{j},:))*Xcol{j};
            end
        end
        % P update
        if P_updt
            for j = 1 : q
                xHat = U(rowIdx{j},:)*V(:,j);
                err = err + norm(Xcol{j} - xHat);
                start = 1;
                for s = 1 : length(r_All{j})
                    stop = start + r_All{j}(s) - 1;           
                    [~, idx1] = sort(xHat(start:stop));
                    [~, idx2] = sort(Xcol{j}(start:stop));
                    idx1 = start - 1 + idx1;
                    idx2 = start - 1 + idx2;
                    rowIdx{j}(idx2) = rowIdx{j}(idx1);
                    start = start + r_All{j}(s);
                end
            end
        end
        % U update
        for t_ = 1 : T_in
            gradU = 0*gradU;
            for j = 1 : q
                gradU(rowIdx{j},:) = gradU(rowIdx{j},:) + (U(rowIdx{j},:)*V(:,j) - Xcol{j})*V(:,j)';
            end
            U = U - eta*gradU;
        end
        %Xhat = U*V;
        objVals(i+1) = err;%norm(Xzeros(idx) - Xhat(idx))/norm(Xzeros(idx));
        if real == 0
            [Uproj,~,~] = qr(U,'econ');
            SDVals(i + 1) = norm(Ustr - Uproj*(Uproj'*Ustr));  
            if mod(i, 5) == 0
                disp(['P_updt = ', num2str(P_updt), '. Iteration = ', num2str(i), ...
                    '. SD = ', num2str(norm(Ustr - Uproj*(Uproj'*Ustr))), ...
                    '. rel Err x = ', num2str(norm(U*V - Ustr*Bstr,'fro')/norm(Ustr*Bstr,'fro')),....
                    ' objVals = ', num2str(err)])
            end
        end 
    end


