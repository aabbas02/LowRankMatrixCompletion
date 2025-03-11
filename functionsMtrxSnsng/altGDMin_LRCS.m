function [SDVals] = altGDMin_MtrxSensing(Ak_, yk_, U0,r,T,Ustr)
    m = size(Ak_{1},1);
    n = size(Ak_{1},2);
    SDVals = zeros(T+1,1);
    U = U0;
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    q = length(yk_);
    B = zeros(r,q);
    gradU = zeros(n,r);
   
    for i = 1 : T
        % V update
        % sort the entries of B
        for k = 1 : q
            B(:,k) = (Ak_{k}*U)\yk_{k};
        end       
        % U update
        X = U*B;
        if i == 1
            X0 = X;
        end
        gradU = 0*gradU;
        for k = 1 : q
            gradU = gradU + Ak_{k}'*(Ak_{k}*X(:,k)-yk_{k})*B(:,k)';
        end
        eta = 0.5/norm(X0)^2;
        U = U - (eta/m)*gradU;
        [U,~,~] = qr(U,'econ');
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
    end
    
end 