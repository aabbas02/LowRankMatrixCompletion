function [SDVals,objVals] = altMinCntrlNew(n,q,r,Uinit, ...
                                            Ustr,T, ...
                                            rowIdx,Xcol,colIdx,Xrow, Xhat0,idx,Xzeros)
    U = Uinit;
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr );
    %[U0Cllps, S0Cllps, V0Cllps] = svds(XzerosCllps,r);
    %U0Cllps = U0Cllps(:,1:r); S0Cllps = S0Cllps(1:r,1:r); V0Cllps = V0Cllps(:,1:r);
    %X0Cllps = U0Cllps*S0Cllps*V0Cllps';
    objVals = zeros(T+1,1);
    objVals(1) = norm(Xzeros(idx) - Xhat0(idx))/norm(Xzeros(idx));   
    V = zeros(r,q);    
    for i = 1 : T
        % V update
        parfor j = 1 : q
            V(:,j) = U(rowIdx{j},:)\Xcol{j};
        end
        % U update 
        parfor j = 1 : n
			U(j,:) = V(:,colIdx{j})'\Xrow{j};
        end
        %[Uproj,~,~] = qr(U,'econ');
        %SDVals(i + 1) = norm( (eye(n) - Uproj*Uproj')*Ustr);  
        Xhat = U*V;
        objVals(i+1) = norm(Xzeros(idx) - Xhat(idx))/norm(Xzeros(idx));
    end
