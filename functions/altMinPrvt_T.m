function [SDVals,timeArr,timeArrComm]  = altMinPrvt_T(Xzeros,r,p, ...
                                        Ustr,T, ...
                                        rowIdx,Xcol,T_inner,numWrkrs,Tsvd)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % SVD initialization
    tStart = tic;
    [Unew,S] = fedSvd(Xzeros/p,r,Tsvd,numWrkrs);
    mu = min(vecnorm(Unew'))*sqrt(n/r);
    % Initialization of Uzero and Vzero after Projection 
    const1 = mu*sqrt(r/n);
    Unew = Unew .* repmat(const1./sqrt(sum(Unew.^2,2)),1,r); % could replace sqrt(sum(Unew.^2,2)) by vecnorm(Unew');
    % 2nd SVD
    U0_init = orth(Unew);
    U = U0_init(:,1:r);    
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    etaIn = 1/(S(1,1)^2*p);
    tmp = cell(q,1);
    timeArr = zeros(T+1,1);
    timeArrComm = zeros(T+1,1);
    V = zeros(r,q);
    tWrkrs1 = zeros(q,1);
    tWrkrs2 = zeros(q,1);
    for i = 1 : T
        % V update
        tStartFed1 = tic;
        parfor j = 1 : q
            tStartWrkr1 = tic;
            V(:,j) = U(rowIdx{j},:)\Xcol{j};
            tWrkrs1(j) = toc(tStartWrkr1);
        end
        tStopFed1 = toc(tStartFed1);
        tWrkrNodeAvg1 = (sum(tWrkrs1)/q)*(q/numWrkrs);
        tStopFed1 = tStopFed1 - tWrkrNodeAvg1;
        % U update
        tFed2Total = 0;
        for tIn = 1 : T_inner
                % node terms
                tStartFed2 = tic;
                parfor j = 1 : q
                    tStartWrkr2 = tic;
                    tmp{j} =  ( U(rowIdx{j},:)*V(:,j) - Xcol{j} )*V(:,j)';
                    tWrkrs2(j) = toc(tStartWrkr2);
                end
                tStopFed2 = toc(tStartFed2);
                tWrkrNodeAvg2 = (sum(tWrkrs2)/q)*(q/numWrkrs);
                tStopFed2 = tStopFed2 - tWrkrNodeAvg2; 
                tFed2Total = tFed2Total + tStopFed2;
                for j = 1 : q
                    U(rowIdx{j},:) = U(rowIdx{j},:) - etaIn*tmp{j}; 
                end
                
        end
        timeEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i)  +  timeEnd;
        timeArrComm(i + 1) = timeArrComm(i) + tStopFed1 + tFed2Total;
        [Uproj,~,~] = qr(U,'econ');
        SDVals(i + 1) = norm( (eye(n) - Uproj*Uproj')*Ustr ,'fro' );
        tStart = tic;
    end
end