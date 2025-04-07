function [SDVals,objVals] = minCllpsInit(n,q,r,r_,Uinit, ...
                                         Ustr,T, ...
                                         rowIdx,Xcol, Xhat0,idx,Xzeros,real)

% why gradient descent for U on the collapsed objective instead of least-squares?
% because, the collapsed matrix is different across all columns. 
% chain Rule: https://tutorial.math.lamar.edu/classes/calciii/chainrule.aspx
U = Uinit;
SDVals = zeros(T+1,1);
Uproj = qr(U);
SDVals(1) =  norm(Ustr - Uproj*(Uproj'*Ustr),'fro'); 
objVals = zeros(T+1,1);
% Step - size calculation based off of operator norm of Xinit
Xinit = zeros(n,q);
numObserved = 0;
for k =  1 : q
    Ucllpsk = zeros(length(r_{k}),r);
    XcolCllpsk = zeros(length(r_{k}),1);
    %1. form collapsed matrix Ucllpsk and collapsed Xcolk
    start = 1;
    for i = 1 : length(r_{k})
        blkSize = r_{k}(i);
        rowsBlock = rowIdx{k}(start: start + blkSize - 1);            
        Ucllpsk(i,:) = sum(U(rowsBlock,:));
        XcolCllpsk(i) = sum(Xcol{k}(start: start + blkSize - 1));
        start = start  + blkSize;
        % Option A for X init
        numObserved = numObserved + 1;
        Xinit(rowsBlock(1),k) = XcolCllpsk(i);
    end
    % Option B for X init
    %numObserved = length(Xcol{k});
    %2. least - squares update
    %Xinit(rowIdx{k},k) = U(rowIdx{k},:)*(pinv(Ucllpsk)*XcolCllpsk);
    % OR
    %Xinit(rowIdx{k},k) = Xcol{k};     
end
p = numObserved/(n*q);
[~,S] = svds(Xinit/p);
eta = 0.9/(S(1,1)^2*p); 
gradU = zeros(n,r);
for t = 1 : T
    gradU = 0*gradU;
    err = 0;
    for k =  1 : q
        Ucllpsk = zeros(length(r_{k}),r);
        XcolCllpsk = zeros(length(r_{k}),1);
        %1. form collapsed matrix Ucllpsk and collapsed Xcolk
        start = 1;
        for i = 1 : length(r_{k})
            blkSize = r_{k}(i);
            rowsBlock = rowIdx{k}(start: start + blkSize - 1);            
            Ucllpsk(i,:) = sum(U(rowsBlock,:));
            XcolCllpsk(i) = sum(Xcol{k}(start: start + blkSize - 1));
            start = start  + blkSize;
        end
        %2. least - squares update
        bk = pinv(Ucllpsk)*XcolCllpsk;
        %3. get gradient with respect to Ucllpsk 
        diffk = Ucllpsk*bk - XcolCllpsk; % vector of cardinality |Omega_k|
        err = err + norm(diffk);
        gradk = diffk*bk'; % this is an |Omega_k| x r matrix
        %4. get gradient with respect to  U -  chain rule
        start = 1;
        for i = 1 : length(r_{k})
            blkSize = r_{k}(i);
            rowsBlock = rowIdx{k}(start: start + blkSize - 1);
            gradU(rowsBlock,:) = gradU(rowsBlock,:) + gradk(i,:); % All rows that make up block i, denoted rowsBlock here,
                                                                  % have the same gradient
            start = start + blkSize;
        end
    end
    U = U - eta*gradU;
    [Uproj,~,~] = qr(U,'econ');
    SDVals(t+1) = norm(Ustr - Uproj*(Uproj'*Ustr));    
    if mod(t,50) == 0 && real == 0
        disp(['ObjVal = ', num2str(err), '. Iter =  ', num2str(t), '. Gradient Norm = ', num2str(norm(gradU)), '. Subspace Distance = ', num2str(SDVals(t+1))])
        if norm(gradU) < 1e-8
            break
        end
    end
end
