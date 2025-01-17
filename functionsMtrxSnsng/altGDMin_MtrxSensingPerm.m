function [SDVals] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_,U0Cllps,r,T,Ustr)
    m = size(Ak_{1},1);
    n = size(Ak_{1},2);
    SDVals = zeros(T+1,1);
    U = U0Cllps;
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    q = length(ykPerm_);
    B = zeros(r,q);
    gradU = zeros(n,r);
    for i = 1 : T
        % V update
        % sort the entries of B
        for k = 1 : q
            % min over b_k
            B(:,k) = (Ak_{k}*U)\ykPerm_{k};
            % min over P_k
            b = ykPerm_{k}; % no need to sort this at every iteration, can possibly do this once
            % Permutation update
            [~,idx1] = sort(B(:,k));
            [~,idx2] = sort(b);
            Ak_{k}(idx2) = Ak_{k}(idx1);
            %b(idx1) = ykPerm_{k}(idx2); 
        end

       
        % U update
        X = U*B;
        if i == 1
            X0 = X;
        end
        
        gradU = 0*gradU;
        for k = 1 : q
            gradU = gradU + Ak_{k}'*(Ak_{k}*X(:,k)-ykPerm_{k})*B(:,k)';
        end
        eta = 0.5/norm(X0)^2;
        U = U - (eta/m)*gradU;
        [U,~,~] = qr(U,'econ');
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
    end
end 