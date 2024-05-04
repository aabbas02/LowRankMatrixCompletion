function [SDVals, timeArr] = altMinCntrl(Xzeros,r,p, ...
                                        Ustr,T, ...
                                        rowIdx,Xcol,colIdx,Xrow,numWrkrs,Tsvd)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % SVD initialization
    tStart = tic;
    [U,~,~] = svds(Xzeros/p,r);
    U = U(:,1:r);
    %[U,~] = fedSvd(Xzeros/p,r,Tsvd,numWrkrs);
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    timeArr = zeros(T+1,1);
    V = zeros(r,q);
    for i = 1 : T
        %tStart = tic;
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