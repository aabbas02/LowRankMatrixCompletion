function [SDVals,objVals] = altGDMin_LRMC(n,q,r,r_,p,Uinit, ...
                                        Ustr,T, ...
                                        rowIdx,Xcol,colIdx,Xrow, Xhat0,idx,Xzeros,real,B_init)
    % Given U^{(0)}, B^{(0)}
    % Update P_j for all j in [q]. This is to to update rowIdx{j}
    % Update b_j for all j in [q]. This is a least-squares problem with
    % P_jU as the sensing matrix.
    % Update U by gradient descent - Why is this easier than updating U by
    % least-squares?
    V = B_init;
    U = Uinit;
    gradU = zeros(n,r);
    [~,S] = svds(Xzeros/p);
    % Initialization of Uzero and Vzero after Projection 
    eta = 0.1/(S(1,1)^2*p); 
    
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
            %[~,idx1] = sort(xHat);
            %[~,idx2] = sort(Xcol{j}); % make sure this does not do an in-place sort
            %rowIdx{j}(idx2) = rowIdx{j}(idx1);
            for s = 1 : length(r_{j})
                start = sum(r_{j}(1:s)) - r_{j}(s) + 1;
                stop = sum(r_{j}(1:s));           
                [~,idx1] = sort(xHat(start:stop));
                [~,idx2] = sort(Xcol{j}(start:stop));
                idx1 = start - 1 + idx1;
                idx2 = start - 1 + idx2;
                rowIdx{j}(idx2) = rowIdx{j}(idx1);
                %Ak_{k}(idx2,:) = Ak_{k}(idx1,:);
            end

        end
        % V update
        for j = 1 : q
            V(:,j) = U(rowIdx{j},:)\Xcol{j};
        end
        % U update
        for T_in = 1:1
            gradU = 0*gradU;
            for j = 1 : q
                gradU(rowIdx{j},:) = gradU(rowIdx{j},:) + (U(rowIdx{j},:)*V(:,j) - Xcol{j})*V(:,j)';
            end
            U = U - eta*gradU;
        end
        if real == 0
            [Uproj,~,~] = qr(U,'econ');
            SDVals(i + 1) = norm(Ustr - Uproj*(Uproj'*Ustr));  
            if mod(i,5) == 0
                norm(Ustr - Uproj*(Uproj'*Ustr))
            end
            %U = Uproj;
        end 
        Xhat = U*V;
        objVals(i+1) = norm(Xzeros(idx) - Xhat(idx))/norm(Xzeros(idx));
    end
