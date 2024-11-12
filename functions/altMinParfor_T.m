function [SDVals, timeArr,timeArrComm] = altMinParfor_T(Xzeros,r,p, ...
                                          Ustr,T, ...
                                          rowIdx,Xcol,colIdx,Xrow,numWrkrs,Tsvd)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % SVD initialization
    tStart = tic;
    [Unew,~] = fedSvd(Xzeros/p,r,Tsvd,numWrkrs);
    mu = min(vecnorm(Unew'))*sqrt(n/r);
    % Initialization of Uzero and Vzero after Projection 
    const1 = mu*sqrt(r/n);
    Unew = Unew .* repmat(const1./sqrt(sum(Unew.^2,2)),1,r); % could replace sqrt(sum(Unew.^2,2)) by vecnorm(Unew');
    % 2nd SVD
    [U0_init,~,~] = svds(Unew,r);
    U = U0_init(:,1:r);
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    timeArr = zeros(T+1,1);
    timeArrComm = zeros(T+1,1);
    V = zeros(r,q);
    tWrkrs = zeros(q,1);
    for i = 1 : T
        % V update
        tStartFed = tic;
        parfor j = 1 : q
            tStartWrkr = tic;
            V(:,j) = U(rowIdx{j},:)\Xcol{j};
            tWrkrs(j) = toc(tStartWrkr);
        end
        tEndFed = toc(tStartFed);
        tWrkrNodeAvg = ( sum(tWrkrs)/q ) * (q/numWrkrs);
        tEndFed = tEndFed - tWrkrNodeAvg;
        % U update 
        for j = 1 : n
			U(j,:) = V(:,colIdx{j})'\Xrow{j};
        end
        timeEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i)  +  timeEnd;
        timeArrComm(i + 1) = timeArrComm(i) + tEndFed;
        [Uproj,~,~] = qr(U,'econ');
        SDVals(i + 1) = norm( (eye(n) - Uproj*Uproj')*Ustr ,'fro' );
        tStart = tic;
    end
end