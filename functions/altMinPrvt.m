function [SDVals,timeArr]  = altMinPrvt(Xzeros,r,p, ...
                                        Ustr,T, ...
                                        rowIdx,Xcol,T_inner,numWrkrs,Tsvd)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % SVD initialization
    %[U,S,~] = svds((1/p)*Xzeros,r);
    %U = U(:,1:r);
    tStart = tic;
    [U,S] = fedSvd(Xzeros/p,r,Tsvd,numWrkrs);
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