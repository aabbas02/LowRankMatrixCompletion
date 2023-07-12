clc
close all
clear all
% n,q,r,p set according to Figure 1a, 1c. 
% For figure 1b, 1d, change r = 10, p = 0.1.
r = 5;
n = 1*3000;
q = 1*5000;
% generate rank-r X*
X = randn(n,q);
[U,Sigma, V] = svd(X,'econ');
VT = V';
X = U(:,1:r)*Sigma(1:r,1:r)*VT(1:r,:);
Ustr = U(:,1:r);
%--------------------------------------------------------------------------
% sub sample X with probaility p
p = 0.05;
idx = randperm(n*q);
idx = idx(1:round(p*n*q));
idxC = setdiff(1:n*q,idx);
[row,col] = ind2sub([n,q],idx);
% Instantiate Y = X_Omega
Xzeros = zeros(n,q); 
Xzeros(idx) = X(idx);
% Mke cell arrays for row and column indices
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
% Gradient descent iterations and Plotting at 1, 1 + space, 1 + 2*space,...
space = 5;
T = 100 + space;
% ------ Start parallelpool (this takes a few minutes)
parpool
%--- AltGD (Benchmark) 
krbstPCA = 2; %number of extra  (compared to AltMin) iterations for AltGDMin
[SDValsrbstPCA, timeArrRbstPCA] = AltGD(Xzeros,r,Ustr, krbstPCA*T,p, idxC);
0
% --- AltGDMin (Proposed)
kAltGDMin = 2; % number of extra iterations for AltGDMin
[rowSort, colSort, ~] = find(Xzeros);
idx = sub2ind([n,q],rowSort,colSort);
Xvec = X(idx);
[SDValsAltGDMin, timeArrAltGDMin] = altGDMin(Xzeros,r,p, ...
                                             Ustr, kAltGDMin*T, ...
                                             rowIdx,Xcol);
1
% --- AltMin (Benchmark)
[SDValsAltMin, timeArrAltMin] = altMin(Xzeros,r,p, ...
									   Ustr,T, ...
									   rowIdx,Xcol,colIdx,Xrow);
2
% --- AltMinPrvt (Benchmark)
T_inner = 10; % number of gradient descent iterations
[SDValsAltMinPrvt, timeArrAltMinPrvt] = altMinPrvt(Xzeros,r,p, ...
                                                   Ustr,T,...
                                                   rowIdx,Xcol,T_inner);
%--------------------------------------------------------------------------
3
% --- ProjGD 
kProjGD = 1;
[SDValsProjGD,timeArrProjGD] = ProjGD(Xzeros,r,Ustr, kProjGD*T,p,idxC);
% --- Plotting (Error against iteration) 
% --- 
t = 1:space:T;
figure
semilogy(t,SDValsAltMin(t),'DisplayName', ...
     'AltMin (Cntrl.)', ...
     'LineWidth',2.25,'Marker','o','MarkerSize',9)
hold on
semilogy(t,SDValsAltGDMin(t),'DisplayName', ...
     'AltGDMin', ...
     'LineWidth',2.25,'Marker','+','MarkerSize',9)
semilogy(t,SDValsrbstPCA(t),'DisplayName', ...
     'AltGD', ...
     'LineWidth',2.25,'Marker','*','MarkerSize',9)
semilogy(t,SDValsProjGD(t),'DisplayName', ...
     'ProjGD', ...
     'LineWidth',2.25,'Marker','diamond','MarkerSize',9)
semilogy(t,SDValsAltMinPrvt(t),'DisplayName',...
		'AltMin (Prvt.)',...
        'LineWidth',1.75,'Marker','square','MarkerSize',9,...
        'Color',"#A2142F")
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
stringTitle = ['ErrAgnstIter_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
               'T_',num2str(T),'_id',num2str(randi(1e3,1)),'.fig'];
 
% save matlab (.fig) figure
savefig(stringTitle)
stringTitlePdf = ['ErrAgnstIter_','_n_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                'T_',num2str(T),'_id',num2str(randi(1e3,1)),'.pdf'];
% save pdf figure
exportgraphics(gcf,stringTitlePdf)
% Plotting (Error against time) %1(b), 1(d)%
figure
t = 1:space:T;
semilogy(timeArrAltMin(t),SDValsAltMin(t), ...
        'DisplayName', 'AltMin (Cntrl.)', ...
         'LineWidth',2.25,'Marker','o','MarkerSize',9)
hold on
semilogy(timeArrAltGDMin(1:space:T),SDValsAltGDMin(1:space:T), ...
        'DisplayName','AltGDMin', ...
        'LineWidth',2.25,'Marker','+','MarkerSize',9)
semilogy(timeArrRbstPCA(1:space:T),SDValsrbstPCA(1:space:T),...
        'DisplayName','AltGD', ...
        'LineWidth',2.25,'Marker','*','MarkerSize',9)
semilogy(timeArrProjGD(1:space:T),SDValsProjGD(1:space:T), ...
        'DisplayName', 'ProjGD', ...
        'LineWidth',2.25,'Marker','diamond','MarkerSize',9)
semilogy(timeArrAltMinPrvt(1:space:T),SDValsAltMinPrvt(1:space:T), ...
        'DisplayName','AltMin (Prvt.)',...
        'LineWidth',2.25,'Marker','square','MarkerSize',9,...
        'Color',"#A2142F")
grid on
% set x limit
idxMax = find(SDValsAltGDMin < 2e-15,1);
tMax = timeArrAltGDMin(idxMax);
xlimits = xlim();
xlimits(2) = min(xlimits(2),20);
xlim(xlimits);
legend('Interpreter','Latex','Fontsize',9)
ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$','Interpreter','Latex','Fontsize',15)
xlabel('Time/seconds',FontSize=11)
title("n = " + n + ", q = " + q +...
      ", r = " + r + ", p = " + p + '.',...
       'Interpreter', 'Latex', 'FontSize',14)
stringTitle = ['ErrAgnstTime',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
               'T_',num2str(T),'_id',num2str(randi(1e3,1)),'.fig'];

% save matlab (.fig) figure
savefig(stringTitle)
stringTitlePdf = ['ErrAgnstTime',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                'T_',num2str(T),'_id',num2str(randi(1e3,1)),'.pdf'];
% save pdf figure
exportgraphics(gcf,stringTitlePdf)
%--- Function routines for Algorithms
function [SDVals, timeArr] = altGDMin(Xzeros,r,p,...
                                      Ustr, T,...
                                      rowIdx,Xcol)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    [U0_init,S,~] = svd((1/p)*Xzeros,'econ');
    U = U0_init(:,1:r);
    eta = 1/(S(1,1)^2*p); %S(1,1)^2 = norm(Xzeros,2)^2/p^2, %eta = p/norm(Xzeros,2)^2;
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    timeArr = zeros(T+1,1);
    grad = cell(q,1);
    gradMat = zeros(n,r);
    for i = 1 : T
        tStart = tic;
        parfor j = 1 : q
            B(:,j) = U(rowIdx{j},:)\Xcol{j};
        end        
        parfor j = 1 : q % could avoid the second parfor loop
            grad{j} = ( U(rowIdx{j},:)*B(:,j) - Xcol{j} )*B(:,j)';
        end
        tCenter = tic;
        gradMat = 0*gradMat;
        for j = 1 : q               
            gradMat(rowIdx{j},:) =  gradMat(rowIdx{j},:) + grad{j};
        end
        U = U - eta*gradMat;
        timeExtra = toc(tCenter);
        % Project U by using SVD or QR
        [U,~] = qr(U,'econ');
        tEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i) + tEnd - timeExtra;
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
    end
end
function [SDVals,timeArr]  = altMinPrvt(Xzeros,r,p, ...
                                        Ustr,T, ...
                                        rowIdx,Xcol,T_inner)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % SVD initialization
    [U,S,~] = svd((1/p)*Xzeros,'econ');
    U = U(:,1:r);
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    etaIn = 1/(S(1,1)^2*p);
    tmp = cell(q,1);
    timeArr = zeros(T+1,1);
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
        tQR = tic;
        [Uproj,~,~] = qr(U,'econ');
        timeQR = toc(tQR);
        timeEnd = toc(tStart);
        %timeExtra
        timeArr(i + 1) = timeArr(i)  +  timeEnd  - timeQR - timeExtra;
        SDVals(i + 1) = norm( (eye(n) - Uproj*Uproj')*Ustr ,'fro' );
    end
end
function [SDVals, timeArr] = altMin(Xzeros,r,p, ...
                                    Ustr,T, ...
                                    rowIdx,Xcol,colIdx,Xrow)
    n = size(Xzeros,1);
    q = size(Xzeros,2);
    % SVD initialization
    [U,~,~] = svd((1/p)*Xzeros,'econ');
    U = U(:,1:r);
    SDVals = zeros(T+1,1);
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    timeArr = zeros(T+1,1);
    for i = 1 : T
        tStart = tic;
        % V update
        parfor j = 1 : q
            V(:,j) = U(rowIdx{j},:)\Xcol{j};
        end
        % U update 
		parfor j = 1 : n
			U(j,:) = V(:,colIdx{j})'\Xrow{j};
		end
        tQR = tic;
        [Uproj,~,~] = qr(U,'econ');
        timeQR = toc(tQR);
        timeEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i)  +  timeEnd  - timeQR;
        SDVals(i + 1) = norm( (eye(n) - Uproj*Uproj')*Ustr ,'fro' );
    end
end
function [SDVals, timeArr] = AltGD(Xzeros,r,Ustr,T,p,idxC)
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
        [U,S,V] = svd(X,'econ');
        U = U(:,1:r);
        S = S(1:r,1:r);
        V = V(:,1:r);
        X = U*S*V';
        timeEnd = toc(tStart);
        timeArr(i + 1) = timeArr(i) +  timeEnd;
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr,'fro' );
    end
end