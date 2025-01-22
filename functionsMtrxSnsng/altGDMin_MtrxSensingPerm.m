function [SDVals] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_,AkCllps_,ykCllps_,U0Cllps,r,T,Ustr)
    %---
    % Algorithm
    % Init: U^(0), B^(0)
    % Steps: min P, min U, min B
    % if t == 0, update bk by collapsed estimate, else update by full
    % estimate
    % repeat
    m = size(Ak_{1},1);
    n = size(Ak_{1},2);
    SDVals = zeros(T+1,1);
    U = U0Cllps;
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    q = length(ykPerm_);
    B = zeros(r,q);
    gradU = zeros(n,r);
    for i = 1 : T
        % b_k update and P_k update in for k = 1 : q loop below
        for k = 1 : q
            % Least-squares B_k update 
            if i == 1 % collapsed least-squares
                B(:,k) = (AkCllps_{k}*U)\ykCllps_{k};
            else     % full measurements least-squares
                B(:,k) = (Ak_{k}*U)\ykPerm_{k};
            end
            % min over P_k
            yHat_k = Ak_{k}*U*B(:,k);
            [~,idx1] = sort(yHat_k);
            [~,idx2] = sort(ykPerm_{k});
            Ak_{k}(idx2,:) = Ak_{k}(idx1,:);
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
        eta = 5e-1/norm(X0)^2;
        U = U - (eta/m)*gradU;
        [U,~,~] = qr(U,'econ');
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
    end
end 