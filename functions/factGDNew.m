function [SDVals, timeArr] = factGDNew(Xzeros,r,Ustr,T,p,rowIdx,colIdx,Xcol,numWrkrs,Tsvd)
    % rowIdx - cell array of size 1 x q
    % colIdx - cell array of size 1 x n
    % Xcol - cell array of size q x 1
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % Initialization using SVD 
    tStart = tic;
    [U,Sig,V] = fedSvd_UV(Xzeros/p,r,Tsvd,numWrkrs);
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
    rowIdx_ = cell(1,numWrkrs);
    colIdx_ = cell(1,numWrkrs);
    Xcol_ = cell(1,numWrkrs);
    V_ = cell(1,numWrkrs);
    for j = 1 : numWrkrs
        offset = (j-1)*q/numWrkrs;
        V_{j} = V(offset + 1: offset + q/numWrkrs,:);
        rowIdx_{j} = cell2mat(rowIdx(offset + 1 : j*q/numWrkrs)');
        cols = zeros(1,length(rowIdx_{j}));
        strt = 0;
        for k = offset + 1 : j*q/numWrkrs
            lenK = length(rowIdx{k});
            cols(strt + 1 : strt + lenK) = k;
            strt = strt + lenK;
        end
        colIdx_{j} = cols;
        Xcol_{j} = cell2mat(Xcol(offset + 1 : j*q/numWrkrs));
    end
    T1_ = zeros(n,r,numWrkrs);
    T3_ = zeros(r,r,numWrkrs);
    for i = 1 : T
        parfor j = 1 : numWrkrs
            offset = (j-1)*q/numWrkrs;
            diff = sum(U(rowIdx_{j},:).*V_{j}(colIdx_{j} - offset,:), 2) - Xcol_{j};
            temp = sparse(rowIdx_{j},colIdx_{j} - offset,diff,n,q/numWrkrs); %(This is the sparse matrix (X - Y)_{Omega})
            T1_(:,:,j) = temp*V_{j};
            T2 = temp'*U;
            V_{j} = V_{j} - steplength*T2/p %GD update w.r.t to loss function/ First Data Exchange
            T3_(:,:,j) = V_{j}'*V_{j};
        end
        % U update - at the center
        T1 = sum(T1_,3);
        T4 = U'*U; % this is also transmitted to nodes in the second data exchange done for updating B
        T3 = sum(T3_,3);
        Unew = U - steplength*T1/p - steplength/16*U*(T4-T3); %U is n x r, V is q x r
        Unew = Unew .* repmat(min(ones(n,1),const1./sqrt(sum(Unew.^2,2))),1,r);
        diffM = T3 - T4;
        parfor j = 1 : numWrkrs
            V_{j} = V_{j} -...
            (steplength/16) * V_{j} * diffM; % GD update w.r.t to norm balancing term/Second Data Exchange
            V_{j} = V_{j}.*repmat(min(ones(q/numWrkrs,1),const2./sqrt(sum(V_{j}.^2,2))),1,r);
        end
        timeArr(i+1) = timeArr(i) + toc(tStart);
        % QR to get subspace distance - not used in next iteration
        [U,~] = qr(Unew,'econ');
        SDVals(i+1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
        U = Unew;
        tStart = tic;
    end
end