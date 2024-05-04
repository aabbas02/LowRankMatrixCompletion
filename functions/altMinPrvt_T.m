function [SDVals,timeArr]  = altMinPrvt_T(Xzeros,r,p, ...
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
    V = zeros(r,q);
    for i = 1 : T
        % V update
        parfor j = 1 : q
            V(:,j) = U(rowIdx{j},:)\Xcol{j};
        end
        for tIn = 1 : T_inner
                % node terms
                parfor j = 1 : q
                    tmp{j} =  ( U(rowIdx{j},:)*V(:,j) - Xcol{j} )*V(:,j)';
                end
                for j = 1 : q
                    U(rowIdx{j},:) = U(rowIdx{j},:) - etaIn*tmp{j}; 
                end
        end
        timeEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i)  +  timeEnd;
        [Uproj,~,~] = qr(U,'econ');
        SDVals(i + 1) = norm( (eye(n) - Uproj*Uproj')*Ustr ,'fro' );
        tStart = tic;
    end
end