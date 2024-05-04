function [SDVals,timeArr] =  altGDMinCntrl(Xzeros,Xvec, r,p,...
                                            row,col,rowIdx,Xcol,...
                                            T,Ustr,numWrkrs,Tsvd)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    SDVals = zeros(T+1,1);
    timeArr = zeros(T+1,1);
    tStart = tic;
    [U0_init,S,~] = svds(Xzeros/p,r);
    U = U0_init(:,1:r);
    eta = 1/(S(1,1)^2*p); 
    B = zeros(r,q);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    for i = 1 : T
        %tStart = tic;
        % B update
        for j = 1 : q
            B(:,j) = U(rowIdx{j},:)\Xcol{j};
        end
        diff = sum( U(row,:).*B(:,col)',2 ) - Xvec;
        %make sparse matrix with entries equal to diff vector, supported on idx   
        S = sparse(row,col,diff,n,q);
        grad = S*B';
        U = U - eta*grad;
        % Project U by using SVD or QR
        [U,~] = qr(U,'econ');
        tEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i) + tEnd;
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
        tStart = tic;
    end
end