function [SDVals, timeArr] = AltGDFedSparse(Xzeros,r,Ustr,T,p,rowIdx,colIdx,Xcol,numWrkrs,Tsvd)
    % rowIdx - cell array of size 1 x q
    % colIdx - cell array of size 1 x n
    % Xcol - cell array of size q x 1
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % Initialization using SVD 
    %[Ustr_,S_,Vstr_]= svds(Xzeros/p,r);
    %U = U(:,1:r);
    tStart = tic;
    [U,Sig,V] = fedSvd_UV(Xzeros/p,r,Tsvd,numWrkrs);
    %[V,~] = fedSvd(Xzeros'/p,r,Tsvd,numWrkrs);
    %sprintf('AltGDFedSparse. Error in left-singular vectors %e', norm((eye(n) - U*U')*Ustr_,'fro'))
    %sprintf('AltGDFedSparse. Error in right-singular vectors %e', norm((eye(q) - V*V')*Vstr_,'fro'))
    %sprintf('Error in singular values %e \n', norm(diag(Sig) - diag(S_),'fro'))
    %norm((eye(q) - V*V')*Vstr,'fro')
    %U = Ustr_;
    %V = Vstr_;
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
    for j = 1 : numWrkrs
        offset = (j-1)*q/numWrkrs;
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
    %gradM_ = zeros(n,q/numWrkrs,numWrkrs); % Total upstream =
    %n*q/numWrkrs*numWrkrs
    % Total upstream of T1_,T_2,T_3 = n*r*numWrks + qr + r^2numWrkrs
    T1_ = zeros(n,r,numWrkrs);
    T2_ = zeros(q/numWrkrs,r,numWrkrs);
    T3_ = zeros(r,r,numWrkrs);
    for i = 1 : T
        %tStart = tic;
        parfor j = 1 : numWrkrs
            offset = (j-1)*q/numWrkrs;
            diff = sum(U(rowIdx_{j},:).*V(colIdx_{j},:), 2) - Xcol_{j};
            %linIdx = sub2ind([n,q/numWrkrs],rowIdx_{j},colIdx_{j} - offset);
            %temp = zeros(n,q/numWrkrs);
            %temp(linIdx) = diff; % n x q/numWrkrs
            temp = sparse(rowIdx_{j},colIdx_{j} - offset,diff,n,q/numWrkrs);
            T1_(:,:,j) = temp*V(offset + 1: offset + q/numWrkrs,:); % n x r matrix
            T2_(:,:,j) = temp'*U;  % q/numWrkrs x r 
            T3_(:,:,j) = V(offset + 1: offset + q/numWrkrs,:)'*V(offset + 1: offset + q/numWrkrs,:);
        end
        T1 = sum(T1_,3);
        T2 = reshape( permute(T2_,[ 1 3 2]) , q , []);
        T3 = sum(T3_,3);
        Unew = U - steplength*T1/p - steplength/16*U*(U'*U-T3); %U is n x r, V is q x r
        Vnew = V - steplength*T2/p - steplength/16*V*(T3-U'*U); 
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