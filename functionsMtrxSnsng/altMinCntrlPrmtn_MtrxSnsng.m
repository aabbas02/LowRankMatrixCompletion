function [SDVals, timeArr] = altMinCntrlPrmtn_MtrxSnsng(Xzeros,r,p, ...
                                        Ustr,T, ...
                                        rowIdx,Xcol,colIdx,Xrow,numWrkrs,Tsvd)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % SVD initialization
    tStart = tic;
    [Unew,~,~] = svds(Xzeros/p,r);
    mu = min(vecnorm(Unew'))*sqrt(n/r);
    % Initialization of Uzero and Vzero after Projection 
    const1 = mu*sqrt(r/n);
    Unew = Unew .* repmat(const1./sqrt(sum(Unew.^2,2)),1,r); % could replace sqrt(sum(Unew.^2,2)) by vecnorm(Unew');
    % 2nd SVD
    U0_init = orth(Unew);
    U = U0_init(:,1:r);
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    timeArr = zeros(T+1,1);
    V = zeros(r,q);
    for i = 1 : T
        % V update
        for j = 1 : q
            V(:,j) = U(rowIdx{j},:)\Xcol{j};
        end
        % U update 
        for j = 1 : n
			U(j,:) = V(:,colIdx{j})'\Xrow{j};
        end
        timeEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i)  +  timeEnd;
        [Uproj,~,~] = qr(U,'econ');
        SDVals(i + 1) = norm( (eye(n) - Uproj*Uproj')*Ustr ,'fro' );
        tStart = tic;
    end
end