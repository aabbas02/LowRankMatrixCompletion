function [SDVals,timeArr] = altGDMin(r,eta_c,...
                                     Ustr,Xzeros, T,p, ...
                                     rowIdx,Xcol,numWrkrs,Tsvd)

    % 4 is the number of nodes
    
    % rows_ is a cell array of size 4 x 1, containing indices of observed
    % entries of columns 1...q/4, q/4+1...q/2, q/2+1...3q/4, 3q/4+1...q.
    % Make rows_ a cell array of size 4 x q/4.

    % cols_ is a cell array of size 4x1 such that cols_{j} \in k, where 
    % k = (j-1)*q/4+1:j*q/4, j \in {1,2,3,4}.
    
    % rowIdx is cell of size qx1,containing the indices of the observed rows
    % of each of the q columns

    n = size(Xzeros,1);
    q = size(Xzeros,2);
    tStart = tic;
    [U,S] = fedSvd(Xzeros/p,r,Tsvd,numWrkrs);
    %[U0_init,S,~] = svds(Xzeros/p,r);
    %U = U0_init(:,1:r);
    eta = eta_c/(S(1,1)^2*p); 
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    timeArr = zeros(T+1,1);
    Xcol_ = cell(numWrkrs,q/numWrkrs);
    rowsJ = reshape(rowIdx,[q/numWrkrs,numWrkrs]);
    rowsJ = rowsJ';
    for j = 1 : numWrkrs
        offset = (j-1)*q/numWrkrs;
        for k = 1 : q/numWrkrs
            Xcol_{j,k} = Xcol{offset+k};
        end
    end
    grad_ = zeros(n,r,numWrkrs);  
    for i = 1 : T
        %tStart = tic;
        %start =  ticBytes(gcp);
        parfor j = 1 : numWrkrs
            tmp = zeros(r,q/numWrkrs);
            %M = zeros(n,q/4); 
            %diff = zeros(n*q,1);
            %rows = zeros(n*q,1);
            %cols = zeros(n*q,1);
            diff = [];
            rows = [];
            cols = [];
            rowIdxj = rowsJ(j,:); 
                                % rowIdxj should be cell aray of size  1 x q/4.
                                % alternatively, rows_{j} can be a cell
                                % array of size 4 x 1,
            XcolJ = Xcol_(j,:);
            strt = 0;
            for k = 1 : q/numWrkrs
               rowIdx_jk = rowIdxj{k};
               numRows = length(rowIdx_jk);
               tmp(:,k) = U(rowIdx_jk,:)\XcolJ{k}
               diff(strt + 1 : strt + numRows) = U(rowIdx_jk,:)*tmp(:,k) - XcolJ{k}; 
               rows(strt + 1 : strt + numRows) = rowIdx_jk;
               cols(strt + 1 : strt + numRows) = k;
               strt = strt + numRows;
            end
            diff = diff(1:strt);
            rows = rows(1:strt);
            cols = cols(1:strt);
            %linIdx = sub2ind([n,q/4],rows,cols);
            %M(linIdx) = diff;
            M = sparse(rows,cols,diff,n,q/numWrkrs);
            grad_(:,:,j) = M*tmp'; 
        end
        %tocBytes(gcp,start)
        grad = sum(grad_,3);
        U = U - eta*grad;
        % Project U by using SVD or QR
        [U,~] = qr(U,'econ');
        tEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i) + tEnd;
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
        tStart = tic;
    end        
end
