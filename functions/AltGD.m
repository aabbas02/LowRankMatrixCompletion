function [SDVals, timeArr] = AltGD(Xzeros,r,Ustr,T,p,idxC,numWrkrs,Tsvd)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % Initialization using SVD 
    tStart = tic;
    %[U,Sig] = fedSvd(Xzeros/p,r,Tsvd,numWrkrs);
    %[V,~] =  fedSvd(Xzeros'/p,r,Tsvd,numWrkrs);
    [U,Sig,V]= svds(Xzeros/p,r);
    %U = U(:,1:r);
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    U = U(:,1:r)*sqrt(Sig(1:r,1:r));
    V = V(:,1:r)*sqrt(Sig(1:r,1:r));
    % Step size selection
    step_const = 0.75;
    steplength = step_const / Sig(1,1);
    % Incoherence parameter mu  
    norm_U = max(vecnorm(U'))*sqrt(n/r);
    norm_B = max(vecnorm(V'))*sqrt(q/r);
    mu = max(norm_U, norm_B);
    % Initialization of Uzero and Vzero after Projection 
    const1 = sqrt(4*mu*r/n)*Sig(1,1);
    const2 = sqrt(4*mu*r/q)*Sig(1,1);
    U = U .* repmat(min(ones(n,1),const1./sqrt(sum(U.^2,2))),1,r);
    V = V .* repmat(min(ones(q,1),const2./sqrt(sum(V.^2,2))),1,r);
    timeArr = zeros(T+1, 1);
    for i = 1 : T
        % Gradient Descent for U and V
        %tStart = tic;
        gradM = U*V' - Xzeros;
        gradM(idxC) = 0;
        Unew = U - steplength * (gradM*V)/p - steplength/16*U*(U'*U-V'*V);
        Vnew = V - steplength * (gradM'*U)/p - steplength/16*V*(V'*V-U'*U);
        % Projection
        Unew = Unew .* repmat(min(ones(n,1),const1./sqrt(sum(Unew.^2,2))),1,r);
        Vnew = Vnew .* repmat(min(ones(q,1),const2./sqrt(sum(Vnew.^2,2))),1,r);
        timeArr(i+1) = timeArr(i) + toc(tStart);
        % QR to get subspace distance - not used in next iteration
        [U,~] = qr(Unew,'econ');
        SDVals(i+1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
        U = Unew;
        V = Vnew;
        tStart = tic;
    end
end