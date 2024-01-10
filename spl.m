clc
close all
%clear all
r = 10;
n = 5000;
q = 10000;
% sub sample X with probaility p
eta_c = 0.75;
p = 0.05;
numWorkers = 4;
space = 25;
T = 50 + space;
delete(gcp('nocreate'))
parpool(numWorkers,'IdleTimeout',Inf)
pool = gcp;
wrkrs = pool.NumWorkers;
MC = 1;
tic 
%{
%1
SDValsAltMinParfor = zeros(MC,T+1);
timeAltMinParfor = zeros(MC,T+1);
%2
SDValsAltMinPrvt = zeros(MC,T+1);
timeAltMinPrvt = zeros(MC,T+1);
%3
SDValsAltMinCntrl = zeros(MC,T+1);
timeAltMinCntrl = zeros(MC, T+1);
%4
kAltGD = 2;
SDValsAltGD = zeros(MC,kAltGD*T+1);
timeAltGD = zeros(MC,kAltGD*T+1);
%5
kProjGD = 2;
SDValsProjGD = zeros(MC,kProjGD*T+1);
timeProjGD = zeros(MC,kProjGD*T+1);
%6
kAltGDMin = 2; % number of extra iterations for AltGDMin
SDValsAltGDMin = zeros(MC,kAltGDMin*T+1);
timeAltGDMin = zeros(MC,kAltGDMin*T+1);
%7
SDValsAltGDMineta1 = zeros(MC,kAltGDMin*T+1);
timeAltGDMineta1 = zeros(MC,kAltGDMin*T+1);
%}
%--------------------------------------------------------------------------
% generate rank-r X*
%U = orth(randn(n,r));
%Bstar = randn(r,q);
%X  = U*Bstar;
%Ustr = U(:,1:r);
for mc = 1 : MC
    U = orth(randn(n,r));
    Bstar = randn(r,q);
    X  = U*Bstar;
    Ustr = U(:,1:r);
    idx = randperm(n*q);
    idx = idx(1:round(p*n*q));
    idxC = setdiff(1:n*q,idx);
    [row,col] = ind2sub([n,q],idx); 
    %Instantiate Y = X_Omega
    Xzeros = zeros(n,q); 
    Xzeros(idx) = X(idx);
    %Xvec = X(idx)';
    % Make cell arrays for row and column indices
    rowIdx = cell(q,1); colIdx = cell(n,1);
    Xcol = cell(q,1); Xrow = cell(n,1);
    parfor j = 1 : q
       rowIdx{j} = row(col==j);
       Xcol{j} =  X(rowIdx{j},j);
    end
    parfor i = 1 : n
       colIdx{i} = col(row==i);
       Xrow{i} = X(i,colIdx{i})';
    end
    % --- AltGDMin(Cntrl.)
    [rowSort, colSort, ~] = find(Xzeros);
    idx = sub2ind([n,q],rowSort,colSort);
    Xvec = X(idx);
    [SDAltGDMinCntrl(mc,:),timeAltGDMinCntrl(mc,:)] = altGDMinCntrl(Xzeros,Xvec,r,p,...
                                                                    rowSort,colSort,rowIdx,Xcol,...
                                                                    T,Ustr);
    %-----------------------------------------------------------------------------------
    
    % --- AltMin (Parfor)
    [SDAltMinParfor(mc,:), timeAltMinParfor(mc,:)] = altMinParfor(Xzeros,r,p, ...
									                              Ustr,T, ...
									                              rowIdx,Xcol,colIdx,Xrow);
    % --- AltMin (Cntrl.)
    [SDAltMinCntrl(mc,:), timeAltMinCntrl(mc,:)] = altMinCntrl(Xzeros,r,p, ...
                                                               Ustr,T, ...
                                                               rowIdx,Xcol,colIdx,Xrow);
    % --- AltMin (Fed./Prvt.)
    T_inner =  10;
    [SDAltMinPrvt(mc,:), timeAltMinPrvt(mc,:)]  = altMinPrvt(Xzeros,r,p, ...
                                                                Ustr,T, ...
                                                                rowIdx,Xcol,T_inner);
    % --- AltGD
    kAltGD = 2;
    [SDAltGD(mc,:), timeAltGD(mc,:)] = AltGD(Xzeros,r,Ustr,kAltGD*T,p,idxC);

    % --- ProjGD
    kProjGD = 2;
    [SDProjGD(mc,:), timeProjGD(mc,:)] = ProjGD(Xzeros,r,Ustr,kProjGD*T,p,idxC);

    % --- AltGDMin
    eta_c = 0.75;
    kAltGDMin = 2;
    [SDAltGDMin(mc,:),timeAltGDMin(mc,:)] = altGDMin(r,eta_c, ...
                                                     Ustr,Xzeros, kAltGDMin*T,p, ...
                                                     rowIdx,Xcol);
    eta_c = 1.00;
    [SDAltGDMineta1(mc,:),timeAltGDMineta1(mc,:)] = altGDMin(r,eta_c, ...
                                                             Ustr,Xzeros, kAltGDMin*T,p, ...
                                                             rowIdx,Xcol);
    mc
end
toc
% 1
SDAltMinParfor = sum(SDAltMinParfor,1)/MC;
timeAltMinParfor = sum(timeAltMinParfor,1)/MC;
% 2
SDAltMinCntrl = sum(SDAltMinCntrl,1)/MC;
timeAltMinCntrl = sum(timeAltMinCntrl,1)/MC;
% 3
SDAltMinPrvt = sum(SDAltMinPrvt,1)/MC;
timeAltMinPrvt = sum(timeAltMinPrvt,1)/MC;
% 4
SDAltGD = sum(SDAltGD,1)/MC;
timeAltGD = sum(timeAltGD,1)/MC;
% 5
SDProjGD = sum(SDProjGD,1)/MC;
timeProjGD = sum(timeProjGD,1)/MC;
% 6
SDAltGDMin = sum(SDAltGDMin,1)/MC;
timeAltGDMin = sum(timeAltGDMin,1)/MC;
% 7
SDAltGDMineta1 = sum(SDAltGDMineta1,1)/MC;
timeAltGDMineta1 = sum(timeAltGDMineta1,1)/MC;
% 8
SDAltGDMinCntrl = sum(SDAltGDMinCntrl,1)/MC;
timeAltGDMinCntrl = sum(timeAltGDMinCntrl,1)/MC;
% --- Plot Subspace Distance against Iteration Figure
plotErrvsIter(space, T, n, q, r, p,...
               SDAltMinCntrl, SDAltMinParfor, SDAltMinPrvt, SDAltGD,...
               SDProjGD, SDAltGDMin, SDAltGDMineta1)
%--- Plot Subspace Distance against Time Figure
idx1 = find(SDAltMinParfor < 10^-14,1);
tidx1 = timeAltMinParfor(idx1);

idx2 = find(SDAltMinCntrl < 10^-14,1);
tidx2 = timeAltMinCntrl(idx2);

idx3 = find(SDAltMinPrvt < 10^-14,1);
tidx3 = timeAltMinPrvt(idx3);

idx4 =  find(SDAltGD < 10^-14,1);
tidx4 = timeAltGD(idx4);

idx5 = find(SDProjGD < 10^-14,1);
tidx5 =  timeProjGD(idx5);

idx6 = find(SDAltGDMin < 10^-14,1);
tidx6 =  timeAltGDMin(idx6);

idx7 = find(SDAltGDMineta1 < 10^-14,1);
tidx7 =  timeAltGDMineta1(idx7);

idx8 = find(SDAltGDMinCntrl < 10^-14,1);
tidx8 =  timeAltGDMinCntrl(idx8);

timeSD = max([tidx1,tidx2,tidx3,tidx4,tidx5,tidx6,tidx7,tidx8]);

figure
semilogy(timeAltMinCntrl(timeAltMinCntrl <= timeSD),...
        SDAltMinCntrl(timeAltMinCntrl <= timeSD),...
        'DisplayName', ...
        'AltMin(Cntrl.)', ...
        'LineWidth',1.45,'Marker','o','MarkerSize',5)
hold on
semilogy(timeAltMinParfor(timeAltMinParfor <= timeSD ),...
         SDAltMinParfor(timeAltMinParfor <= timeSD),'DisplayName', ...
        'AltMin(Fed./NotPrvt.)', ...
        'LineWidth',1.45,'Marker','o','MarkerSize',5)
semilogy(timeAltMinPrvt(timeAltMinPrvt <= timeSD ),...
         SDAltMinPrvt(timeAltMinPrvt <= timeSD),'DisplayName', ...
        'AltMin(Fed./Prvt.)', ...
        'LineWidth',1.45,'Marker','o','MarkerSize',5)
semilogy(timeAltGD(timeAltGD <= timeSD ),...
         SDAltGD(timeAltGD <= timeSD),'DisplayName', ...
        'AltGD', ...
        'LineWidth',1.45,'Marker','x','MarkerSize',5)
semilogy(timeProjGD(timeProjGD <= timeSD ),...
         SDProjGD(timeProjGD <= timeSD),'DisplayName', ...
        'ProjGD', ...
        'LineWidth',1.45,'Marker','+','MarkerSize',5)
semilogy(timeAltGDMin(timeAltGDMin <= timeSD ),...
         SDAltGDMin(timeAltGDMin <= timeSD ),'DisplayName', ...
        'AltGDMin(Fed.) $\eta = 3/4$', ...
        'LineWidth',1.45,'Marker','square','MarkerSize',5)
semilogy(timeAltGDMineta1(timeAltGDMineta1 <= timeSD ),...
         SDAltGDMineta1(timeAltGDMineta1 <= timeSD ),'DisplayName', ...
        'AltGDMin(Fed.) $\eta = 1.0$', ...
        'LineWidth',1.45,'Marker','square','MarkerSize',5)
semilogy(timeAltGDMinCntrl(timeAltGDMinCntrl <= timeSD ),...
         SDAltGDMinCntrl(timeAltGDMinCntrl <= timeSD ),'DisplayName', ...
        'AltGDMin(Cntrl.)', ...
        'LineWidth',1.45,'Marker','square','MarkerSize',5)
grid on

legend('Interpreter','Latex','Fontsize',9)
ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$','Interpreter','Latex','Fontsize',15)
xlabel('Time/seconds',FontSize=11)
title("n = " + n + ", q = " + q +...
      ", r = " + r + ", p = " + p + '.',...
       'Interpreter', 'Latex', 'FontSize',14)
cores = feature('numCores');
stringTitle = ['Time_',num2str(wrkrs),'Max',num2str(cores),'_eta_',num2str(eta_c),'_MC_',num2str(MC),...
                  '_n_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                   'T_',num2str(T),'_id',num2str(randi(1e3,1))];

savefig([stringTitle,'.fig']);
% save pdf figure
exportgraphics(gcf,[stringTitle,'.pdf'])


%--- Function routines for Algorithms
function [SDVals,timeArr] =  altGDMinCntrl(Xzeros,Xvec, r,p,...
                                            row,col,rowIdx,Xcol,...
                                            T,Ustr)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    SDVals = zeros(T+1,1);
    timeArr = zeros(T+1,1);
    [U0_init,S,~] = svds(Xzeros/p,r);
    U = U0_init(:,1:r);
    eta = 0.75/(S(1,1)^2*p); 
    B = zeros(r,q);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    for i = 1 : T
        tStart = tic;
        % B update
        parfor j = 1 : q
            B(:,j) = U(rowIdx{j},:)\Xcol{j};
            %rowIdx = row(col==j);
            %B(:,j) = U(rowIdx,:)\X(rowIdx,j);
        end
        diff = sum( U(row,:).*B(:,col)',2 ) - Xvec;
        %make sparse matrix with entries equal to diff vector, supported on idx   
        S = sparse(row,col,diff,n,q);
        grad = S*B';
        U = U - eta*grad;
        % Project U by using SVD or QR
        [U,~] = qr(U,'econ');
        tEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i) + tEnd;
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
    end
end

function [SDVals,timeArr] = altGDMin(r,eta_c,...
                                            Ustr,Xzeros, T,p, ...
                                            rowIdx,Xcol)

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
    [U0_init,S,~] = svds(Xzeros/p,r);
    U = U0_init(:,1:r);
    eta = eta_c/(S(1,1)^2*p); 
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    timeArr = zeros(T+1,1);
    Xcol_ = cell(4,q/4);
    rowsJ = reshape(rowIdx,[q/4,4]);
    rowsJ = rowsJ';
    for j = 1 : 4
        offset = (j-1)*q/4;
        for k = 1 : q/4
            Xcol_{j,k} = Xcol{offset+k};
        end
    end
    grad_ = zeros(n,r,4);  
    for i = 1 : T
        tStart = tic;
        %start =  ticBytes(gcp);
        parfor j = 1 : 4
            tmp = zeros(r,q/4);
            %M = zeros(n,q/4); 
            diff = [];
            rows = [];
            cols = [];
            rowIdxj = rowsJ(j,:); 
                                % rowIdxj should be cell aray of size  1 x q/4.
                                % alternatively, rows_{j} can be a cell
                                % array of size 4 x 1,
            XcolJ = Xcol_(j,:);
            strt = 0;
            for k = 1 : q/4
               rowIdx_jk = rowIdxj{k};
               numRows = length(rowIdx_jk);
               tmp(:,k) = U(rowIdx_jk,:)\XcolJ{k}
               diff(strt + 1 : strt + numRows) = U(rowIdx_jk,:)*tmp(:,k) - XcolJ{k}; 
               rows(strt + 1 : strt + numRows) = rowIdx_jk;
               cols(strt + 1 : strt + numRows) = k;
               strt = strt + numRows;
            end
            %linIdx = sub2ind([n,q/4],rows,cols);
            %M(linIdx) = diff;
            M = sparse(rows,cols,diff,n,q/4);
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
    end        
end
function [SDVals, timeArr] = altMinParfor(Xzeros,r,p, ...
                                    Ustr,T, ...
                                    rowIdx,Xcol,colIdx,Xrow)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % SVD initialization
    [U,~,~] = svds((1/p)*Xzeros,r);
    U = U(:,1:r);
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    timeArr = zeros(T+1,1);
    V = zeros(r,q);
    for i = 1 : T
        tStart = tic;
        % V update
        parfor j = 1 : q
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
    end
end
function [SDVals,timeArr]  = altMinPrvt(Xzeros,r,p, ...
                                        Ustr,T, ...
                                        rowIdx,Xcol,T_inner)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % SVD initialization
    [U,S,~] = svds((1/p)*Xzeros,r);
    U = U(:,1:r);
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    etaIn = 1/(S(1,1)^2*p);
    tmp = cell(q,1);
    timeArr = zeros(T+1,1);
    V = zeros(r,q);
    for i = 1 : T
        timeExtra = 0;
        tStart = tic;
        % V update
        parfor j = 1 : q
            V(:,j) = U(rowIdx{j},:)\Xcol{j};
        end
        for tIn = 1 : T_inner
                % node terms
                parfor j = 1 : q
                    tmp{j} =  ( U(rowIdx{j},:)*V(:,j) - Xcol{j} )*V(:,j)';
                end
                tCenter = tic; 
                for j = 1 : q
                    U(rowIdx{j},:) = U(rowIdx{j},:) - etaIn*tmp{j}; 
                end
                timeExtra =  timeExtra + toc(tCenter);
        end
        timeEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i)  +  timeEnd - 0*timeExtra;
        [Uproj,~,~] = qr(U,'econ');
        SDVals(i + 1) = norm( (eye(n) - Uproj*Uproj')*Ustr ,'fro' );
    end
end
function [SDVals, timeArr] = altMinCntrl(Xzeros,r,p, ...
                                        Ustr,T, ...
                                        rowIdx,Xcol,colIdx,Xrow)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % SVD initialization
    [U,~,~] = svds((1/p)*Xzeros,r);
    U = U(:,1:r);
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    timeArr = zeros(T+1,1);
    V = zeros(r,q);
    for i = 1 : T
        tStart = tic;
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
    end
end
function [SDVals, timeArr] = AltGD(Xzeros,r,Ustr,T,p,idxC)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % Initialization using SVD 
    [U,Sig,V]= svds(Xzeros/p,r);
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
    timeArr = zeros(T+1, 1);
    for i = 1 : T
        % Gradient Descent for U and V
        tStart = tic;
        gradM = U*V' - Xzeros;
        gradM(idxC) = 0;
        Unew = U - steplength * (gradM*V)/p - steplength/16*U*(U'*U-V'*V);
        Vnew = V - steplength * (gradM'*U)/p - steplength/16*V*(V'*V-U'*U);
        % Projection
        Unew = Unew .* repmat(min(ones(n,1),const1./sqrt(sum(Unew.^2,2))),1,r);
        Vnew = Vnew .* repmat(min(ones(q,1),const2./sqrt(sum(Vnew.^2,2))),1,r);
        timeArr(i+1) = timeArr(i) + toc(tStart);
        % QR to get subspace distance - not used in next iteration
        [U,~] = qr(Unew,'econ');
        SDVals(i+1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
        U = Unew;
        V = Vnew;
    end
end
function [SDVals, timeArr] = ProjGD(Xzeros,r,Ustr,T,p,idxC)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    SDVals = zeros(T + 1, 1);
    SDVals(1) = norm(eye(n)*Ustr ,'fro');
    % Initialization
    X = zeros(n,q);
    stepSize = 1/p;
    timeArr = zeros(T+1,1);
    for i = 1 : T
        grad = X - Xzeros;
        grad(idxC) = 0;
        X = X - stepSize*grad;
        tStart = tic;
        [U,S,V] = svds(X,r);
        U = U(:,1:r);
        S = S(1:r,1:r);
        V = V(:,1:r);
        X = U*S*V';
        timeEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i) +  timeEnd;
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr,'fro' );
    end
end
    

function [] = plotErrvsIter(space, T, n, q, r, p,...
                      SDValsAltMinCntrl, SDValsAltMinParfor, SDValsAltMinPrvt, SDValsAltGD,...
                      SDValsProjGD, SDValsAltGDMin, SDValsAltGDMineta1)
                     
    t = 1:space:T;
    figure
    semilogy(t,SDValsAltMinCntrl(t),'DisplayName', ...
        'AltMin(Cntrl).', ...
        'LineWidth',2.25,'Marker','+','MarkerSize',9)
    hold on
    semilogy(t,SDValsAltMinParfor(t),'DisplayName', ...
        'AltMin(Fed./NotPrvt.)', ...
        'LineWidth',2.25,'Marker','square','MarkerSize',9)
    semilogy(t,SDValsAltMinPrvt(t),'DisplayName', ...
        'AltMin(Fed./Prvt.)', ...
        'LineWidth',2.25,'Marker','square','MarkerSize',9)
    semilogy(t,SDValsAltGD(t),'DisplayName', ...
        'AltGD', ...
        'LineWidth',2.25,'Marker','square','MarkerSize',9)
    semilogy(t,SDValsProjGD(t),'DisplayName', ...
        'ProjGD', ...
        'LineWidth',2.25,'Marker','square','MarkerSize',9)
    semilogy(t,SDValsAltGDMin(t),'DisplayName', ...
        'AltGDMin(Fed.), $\eta=3/4$', ...
        'LineWidth',2.25,'Marker','square','MarkerSize',9)
    semilogy(t,SDValsAltGDMineta1(t),'DisplayName', ...
        'AltGDMin(Fed.), $\eta = 1$', ...
        'LineWidth',2.25,'Marker','square','MarkerSize',9)
    
    grid on
    labelsX =cell(1,length(t));
    for i = 1 : length(t)
        labelsX{1,i} = num2str(t(i)-1);
    end
    xticks(t)
    xticklabels( labelsX')
    xlim([0 T])
    legend('Interpreter','Latex')
    ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$','Interpreter','Latex','Fontsize',15)
    xlabel('Iterations t',FontSize=11)
    title("n = " + n + ", q = " + q +...
        ", r = " + r + ", p = " + p + '.',...
        'Interpreter', 'Latex', 'FontSize',15)
%cores = feature('numCores');
% -------------------------------------------------------------------------
end