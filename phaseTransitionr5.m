clc
close all
clear all
warning('off','all')
% To do: incoherence specification in low rank matrix generation
r = 5;
n = 3000;
q = 3000;
% generate rank-r X
%_-----`-----------------------------------
MC = 15;
p_AltMin = [0.003:0.001:0.01]; 
%p_AltGDMin = [0.01:0.005:0.02];
%p_Rest = [0.035:0.005:0.060];
p_Rest = [0.01:0.005:0.05];
%p_ = unique([p_AltMin,p_AltGDMin,p_Rest]);
p_ = unique([p_AltMin,p_Rest]);
numPass = zeros(4,length(p_));
%--------------------------------------------------------------------------
% sub sample X with probaility p
for l = 1 : length(p_)
    p = p_(l);
    tic
    for t = 1 : MC
        X = randn(n,q);
        [U,Sigma, V] = svd(X,'econ');
        VT = V';
        X = U(:,1:r)*Sigma(1:r,1:r)*VT(1:r,:);
        Ustr = U(:,1:r);
        idx = randperm(n*q);
        idx = idx(1:round(p*n*q));
        idxC = setdiff(1:n*q,idx);
        [row,col] = ind2sub([n,q],idx);
        % -------------------------------------------------------------------------
        Xzeros = zeros(n,q);
        Xzeros(idx) = X(idx);
        %--- make cell arrays
        rowIdx = cell(q,1); colIdx = cell(n,1);
        Xcol = cell(q,1); Xrow = cell(n,1);
        for j = 1 : q
            rowIdx{j} = row(col==j);
            Xcol{j} =  X(rowIdx{j},j);
        end
        for i = 1 : n
            colIdx{i} = col(row==i);
            Xrow{i} = X(i,colIdx{i})';
        end
        % GD iterations and Plotting
        space = 10;
        T = 200;
        %--- Robust PCA --------------------------------------------------------------
        %SDValsrbstPCA = zeros(1,T+1);
        [SDValsrbstPCA] = robustPCA(Xzeros,r,Ustr, T,p, idxC);
        numPass(1,l) = numPass(1,l) + min(SDValsrbstPCA < 1e-11);
        %0
        % --- AltGDMin (Decentralized)
        [rowSort, colSort, ~] = find(Xzeros);
        idx = sub2ind([n,q],rowSort,colSort);
        Xvec = X(idx);
        [SDValsAltGDMin] = altGDMinSparseQR(Xvec,r, ...
                                            rowSort,colSort,Ustr,Xzeros, T,p, ...
                                            rowIdx,Xcol);
        numPass(2,l) = numPass(2,l) + min(SDValsAltGDMin < 1e-11);
        %SDValsAltGDMin
        % --- AltMinCntrl
        [SDValsAltMin] = altMinCntrl(Xzeros,r,p, ...
                                     Ustr,T, ...
                                     rowIdx,Xcol,colIdx,Xrow);
        numPass(3,l) = numPass(3,l) + min(SDValsAltMin < 1e-11);
        %-------------------------------------------------------------
        %2
        % --- ProjGD
        [SDValsProjGD] = ProjGD(Xzeros,r,Ustr, T,p,idxC);
        numPass(4,l) = numPass(4,l) + min(SDValsProjGD < 1e-11);
    end
    toc
    l
end
numPass = numPass/MC;
% --- Plotting 
figure
plot(p_,numPass(3,:),'DisplayName', ...
     'AltMin', ...
     'LineWidth',2.25,'Marker','o','MarkerSize',9)
hold on
plot(p_,numPass(2,:),'DisplayName', ...
     'AltGDMin', ...
     'LineWidth',2.25,'Marker','+','MarkerSize',9)
plot(p_,numPass(1,:),'DisplayName', ...
     'AltGD', ...
     'LineWidth',2.25,'Marker','*','MarkerSize',9)
plot(p_,numPass(4,:),'DisplayName', ...
     'ProjGD', ...
     'LineWidth',2.25,'Marker','diamond','MarkerSize',9)
grid on

legend('Interpreter','Latex')
ylabel('$\Pr[\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*) \leq 10^{-10}]$',...
        'Interpreter','Latex','Fontsize',15)
xlabel('Sampling probability $p$','Interpreter','Latex',FontSize=11)
title("n = " + n + ", q = " + q +...
      ", r = " + r +'.',...
       'Interpreter', 'Latex', 'FontSize',15)
stringTitle = ['PhaseTransition_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
               'T_',num2str(T),'_id',num2str(randi(1e3,1)),'.fig'];
savefig(stringTitle)
stringTitlePdf = ['PhaseTransition_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                'T_',num2str(T),'_id',num2str(randi(1e3,1)),'.pdf'];
exportgraphics(gcf,stringTitlePdf)

function [SDVals] = altGDMinSparseQR(Xvec,r, ...
                                     row,col,Ustr,Xzeros, T,p, ...
                                     rowIdx,Xcol)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    [U0_init,S,~] = svd((1/p)*Xzeros,'econ');
    U = U0_init(:,1:r);
    eta = 1/(S(1,1)^2*p); %S(1,1)^2 = norm(Xzeros,2)^2/p^2, %eta = p/norm(Xzeros,2)^2;
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    B = zeros(r,q);
    for i = 1 : T
        for j = 1 : q
            B(:,j) = U(rowIdx{j},:)\Xcol{j};
        end        
        diff = sum(U(row,:).*B(:,col)',2) - Xvec;
        %make sparse matrix with entries equal to diff vector, supported on idx 
        S = sparse(row,col,diff,n,q);
        grad = S*B';
        U = U - eta*grad;
        % Project U by using SVD or QR
        [U,~] = qr(U,'econ');
        SDVals(i+1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
        if SDVals(i+1) < 1e-11
            SDVals = 1e-12;
            break
        end
    end
end
function [SDVals] = altMinCntrl(Xzeros,r,p, ...
                                Ustr,T, ...
                                rowIdx,Xcol,colIdx,Xrow)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    SDVals = zeros(T+1,1);
    % SVD initialization
    [U,~,~] = svd((1/p)*Xzeros,'econ');
    U = U(:,1:r);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
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
        [Uproj,~,~] = qr(U,'econ');
        SDVals(i+1) = norm( (eye(n) - Uproj*Uproj')*Ustr ,'fro' );
        if SDVals(i+1) < 1e-11
            SDVals = 1e-12;
            break
        end
    end
end
function [SDVals] = ProjGD(Xzeros,r,Ustr,T,p,idxC)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    SDVals = zeros(T + 1, 1);
    SDVals(1) = norm(eye(n)*Ustr ,'fro');
    % Initialization
    X = zeros(n,q);
    stepSize = 1/p;
    for i = 1 : T
        grad = X - Xzeros;
        grad(idxC) = 0;
        X = sparse(X - stepSize*grad);
        try
            [U,S,V] = svds(X,r);
        catch
            warning('SVD Error.')
            SDVals = 1e1;
            break
        end
        X = U*S*V';
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr,'fro' );
        if SDVals(i+1) < 1e-11
            SDVals = 1e-12;
            break
        end
    end
end
function [SDVals] = robustPCA(Xzeros,r,Ustr, T,p,idxC)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % Initialization using SVD 
    [U,Sig,V]= svd(Xzeros/p,'econ');
    U = U(:,1:r);
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
    for i = 1 : T
        % Gradient Descent for U and V
        gradM = U*V' - Xzeros;
        gradM(idxC) = 0;
        Unew = U - steplength * (gradM*V)/p - steplength/16*U*(U'*U-V'*V);
        Vnew = V - steplength * (gradM'*U)/p - steplength/16*V*(V'*V-U'*U);
        % Projection
        Unew = Unew .* repmat(min(ones(n,1),const1./sqrt(sum(Unew.^2,2))),1,r);
        Vnew = Vnew .* repmat(min(ones(q,1),const2./sqrt(sum(Vnew.^2,2))),1,r);
        % QR to get subspace distance - not used in next iteration
        [U,~] = qr(Unew,'econ');
        SDVals(i+1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
        if SDVals(i+1) < 1e-11
            SDVals = 1e-12;
            break
        end
        U = Unew;
        V = Vnew;
    end
end
