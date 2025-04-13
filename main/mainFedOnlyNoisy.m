clc
close all
clear all
dir = pwd;
% For linux, replace '\' with '/'
idcs   = strfind(dir,'\');
newdir = dir(1:idcs(end)-1);
cd (newdir)
addpath(genpath('.\functions'));
cd(dir)    
%---------------------------------
r = 20;
n = 5000;
q = 10000;
% sub sample X with probaility p
p = 0.05;
numWrkrs = 10;
space = 25;
T = 25 + space;
delete(gcp('nocreate'))
parpool(numWrkrs,'IdleTimeout',Inf,'SpmdEnabled',false)
pool = gcp;
MC = 25;
tic 
ID = randi(1e3)+5
%ID = 7891359
%ID = 2
%--------------------------------------------------------------------------
% generate rank-r X*
U = orth(randn(n,r));
Ustr = U(:,1:r);
Bstar = randn(r,q);
X  = U*Bstar;
% add noise N
%N = 1e0*randn(n,q)/(sqrt(r)*cond(X)^3);
%noiseVar = 1e-7;
noiseVar = 1/( sqrt(r)*cond(X)^3 );
N = sqrt(noiseVar)*randn(n,q);
X  =  X + N;
Tsvd = 15;
idx = randperm(n*q);
idx = idx(1:round(p*n*q));
idxC = setdiff(1:n*q,idx);
[row,col] = ind2sub([n,q],idx);
%Instantiate Y = X_Omega
Xzeros = zeros(n,q);
Xzeros(idx) = X(idx);
% Make cell arrays for row and column indices
rowIdx = cell(q,1); colIdx = cell(n,1);
Xcol = cell(q,1); Xrow = cell(n,1);
parfor j = 1 : q
    rowIdx{j} = row(col==j);
    Xcol{j} =  X(rowIdx{j},j);
    if j <= n
        colIdx{j} = col(row==j)
        Xrow{j} = X(j,colIdx{j})';
    end
end
saveName = "n_" + num2str(n) + "_q_" + num2str(q) + "_r_" + num2str(r) + "_p_"+...
            num2str(p) + "_MC_" + num2str(MC) + "_randID_" + num2str(ID) + ".mat";
for mc = 1 : MC
    idx = randperm(n*q);
    idx = idx(1:round(p*n*q));
    idxC = setdiff(1:n*q,idx);
    [row,col] = ind2sub([n,q],idx); 
    %Instantiate Y = X_Omega
    Xzeros = zeros(n,q); 
    Xzeros(idx) = X(idx);
    % Make cell arrays for row and column indices
    rowIdx = cell(q,1); colIdx = cell(n,1);
    Xcol = cell(q,1); Xrow = cell(n,1);
    parfor j = 1 : q
       rowIdx{j} = row(col==j);
       Xcol{j} =  X(rowIdx{j},j);
       if j <= n
           colIdx{j} = col(row==j)
           Xrow{j} = X(j,colIdx{j})';
       end
    end
    %-----------------------------------------------------------------------------------
    % --- AltGD (Federated)
    %kAltGD = 2;
    %[SDAltGDFedHalf(mc,:), timeAltGDFedHalf(mc,:)] = factGDEta(Xzeros,r,Ustr,kAltGD*T,p,...
    %                                                 rowIdx,colIdx,Xcol,numWrkrs,Tsvd,0.5);
    % --- AltGD (Federated)
    kAltGD = 2;
    [SDAltGDFed(mc,:), timeAltGDFed(mc,:)] = factGDEta(Xzeros,r,Ustr,kAltGD*T,p,...
                                                       rowIdx,colIdx,Xcol,numWrkrs,Tsvd,0.75);
    % --- AltGD (Federated)
    %kAltGD = 2;
    %[SDAltGDFed1(mc,:), timeAltGDFed1(mc,:)] = factGDEta(Xzeros,r,Ustr,kAltGD*T,p,...
    %                                                     rowIdx,colIdx,Xcol,numWrkrs,Tsvd,1);
    % --- AltMin (Parfor)
    [SDAltMinParfor(mc,:), timeAltMinParfor(mc,:)] = altMinParfor_T(Xzeros,r,p, ...
									                              Ustr,T, ...
									                              rowIdx,Xcol,colIdx,Xrow,numWrkrs,Tsvd);  
    
    % --- AltMin (Fed./Prvt.)
    %T_inner =  10;
    %[SDAltMinPrvt(mc,:), timeAltMinPrvt(mc,:)] = altMinPrvt_T(Xzeros,r,p, ...
    %                                                           Ustr,T, ...
    %                                                           rowIdx,Xcol,T_inner,numWrkrs,Tsvd);

    % --- AltGDMin
    kAltGDMin = 2;
    eta_c = 1.00;
    [SDAltGDMineta1(mc,:),timeAltGDMineta1(mc,:)] = altGDMin_T(r,eta_c, ...
                                                               Ustr,Xzeros, kAltGDMin*T,p, ...
                                                               rowIdx,Xcol,numWrkrs,Tsvd);
    % --- 
    %if  (mod(mc,1) == 0)

   mc
end
toc
save(saveName,...
    "SDAltGDMineta1","timeAltGDMineta1",...
    "SDAltMinParfor","timeAltMinParfor",...
    "SDAltGDFed","timeAltGDFed",...
    "n","p","q","r","numWrkrs","mc","T",...
    "numWrkrs");    %"SDAltMinPrvt","timeAltMinPrvt",...  %"SDAltGDFedHalf","timeAltGDFedHalf",...    "SDAltGDFed1","timeAltGDFed1",...

   

% 1
MC = mc;
SDAltMinParfor = sum(SDAltMinParfor,1)/MC;
timeAltMinParfor = sum(timeAltMinParfor,1)/MC;
% 3
%SDAltMinPrvt = sum(SDAltMinPrvt,1)/MC;
%timeAltMinPrvt = sum(timeAltMinPrvt,1)/MC;
% 4
%SDAltGDFedHalf = sum(SDAltGDFedHalf,1)/MC;
%timeAltGDFedHalf = sum(timeAltGDFedHalf,1)/MC;
% 4
SDAltGDFed = sum(SDAltGDFed,1)/MC;
timeAltGDFed = sum(timeAltGDFed,1)/MC;
%---------------------------------------------
%SDAltGDFed1 = sum(SDAltGDFed1,1)/MC;
%timeAltGDFed1 = sum(timeAltGDFed1,1)/MC;
% 7
SDAltGDMineta1 = sum(SDAltGDMineta1,1)/MC;
timeAltGDMineta1 = sum(timeAltGDMineta1,1)/MC;
%---------------------------------------
pltFedOnly(timeAltMinParfor, SDAltMinParfor, ...
            0, 0, 0, 0, ...
            timeAltGDFed,SDAltGDFed, 0,0,...
            timeAltGDMineta1, SDAltGDMineta1,...
            n, q, r, p, numWrkrs, MC, T,Tsvd,noiseVar)

saveFigFed(1,dir,21,saveName,noiseVar)
%---
function pltFedOnly(timeAltMinParfor, SDAltMinParfor, ...
                    timeAltMinPrvt, SDAltMinPrvt, timeAltGDFedHalf, SDAltGDFedHalf, ...
                    timeAltGDFed,SDAltGDFed, timeAltGDFed1,SDAltGDFed1,...
                    timeAltGDMineta1, SDAltGDMineta1,...
                    n, q, r, p, numWrkrs, MC, T,Tsvd,noiseVar_)

    idx1 = find(SDAltMinParfor < 10^-14,1);
    tidx1 = timeAltMinParfor(idx1);
    %
    %idx3 = find(SDAltMinPrvt < 10^-14,1);
    %tidx3 = timeAltMinPrvt(idx3);
    %
    idx4 = find(SDAltGDFedHalf  < 10^-15,1);
    tidx4 = timeAltGDFedHalf(idx4);
    %
    idx5 = find(SDAltGDFed  < 10^-15,1);
    tidx5 = timeAltGDFed(idx5);
    %
    idx6 = find(SDAltGDFed1  < 10^-15,1);
    tidx6 = timeAltGDFed1(idx6);
    %
    idx7 = find(SDAltGDMineta1 < 10^-14,1);
    tidx7 =  timeAltGDMineta1(idx7);
    %
    timeSD = max([tidx1,tidx5,tidx7]);
    %
    figure;
    semilogy(timeAltMinParfor, ...
        SDAltMinParfor, ...
        'DisplayName', 'AltMin(Fed./NotPrvt.)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 5);
    hold on;
    %semilogy(timeAltMinPrvt(timeAltMinPrvt <= timeSD), ...
    %    SDAltMinPrvt(timeAltMinPrvt <= timeSD), ...
    %    'DisplayName', 'AltMin(Fed./Prvt.)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 5);
    %
    %semilogy(timeAltGDFedHalf(timeAltGDFedHalf <= timeSD), ...
    %    SDAltGDFedHalf(timeAltGDFedHalf <= timeSD), ...
    %    'DisplayName', 'FactGD (Fed. $c = 0.5$)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    %
    semilogy(timeAltGDFed, ...
        SDAltGDFed, ...
        'DisplayName', 'FactGD (Fed. $c = 0.75$)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    %
    %semilogy(timeAltGDFed1(timeAltGDFed1 <= timeSD), ...
    %    SDAltGDFed1(timeAltGDFed1 <= timeSD), ...
    %    'DisplayName', 'FactGD (Fed. $c = 1.00$)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    % 
    semilogy(timeAltGDMineta1, ...
        SDAltGDMineta1, ...
        'DisplayName', 'AltGDMin(Fed.)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 5);
    % 
    grid on;
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$', 'Interpreter', 'Latex', 'Fontsize', 15);
    xlabel('Time/seconds', 'FontSize', 11);
        
    %title("n = " + n + ", q = " + q +...
    %    ", r = " + r + ", p = " + p + '.', ...
    %    'Interpreter', 'Latex', 'FontSize', 14);
    %title("$n = $" +  n + ",$q = $" +  q +...
    %    ",$r = $" +  r + ",$p = $" +  p + "$N \sim \mathcal{N}(0, $" + noiseVar_ + "$)$.", ...
    %    'Interpreter', 'Latex', 'FontSize', 14);
    %
    %cores = feature('numCores');
    %stringTitle = ['FedOnly_', num2str(numWrkrs),'_MC_', num2str(MC), ...
    %    '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_p_', num2str(p), ...
    %    'T_', num2str(T),'_Tsvd',num2str(Tsvd),'_id', num2str(randi(1e3, 1))];
    
    %savefig([stringTitle, '.fig']);
end
%---
function saveFigFed(numData,dir,timeSD,saveName,noiseVar_)
    cd (dir)
    clc
    close all
    %load("n_"+num2str(n)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_p_"+num2str(p),...
    %     +"_MC_"+num2str(mc) + "_randID_" + num2str(ID) + ".mat")
    %load("1.mat");
    load(saveName)
    %T = 0; % You need to define the value of T, as it is used in the code
    numWrkrs = 10; % Define the number of workers
    % --
    SDAltMinParfor_ = zeros(1e3,T+1);
    timeAltMinParfor_ = zeros(1e3,T+1);
    %
    %SDAltMinPrvt_ = zeros(1e3,T+1);
    %timeAltMinPrvt_ = zeros(1e3,T+1);
    %
    %SDAltGDFedHalf_ = zeros(1e3,2*T+1);
    %timeAltGDFedHalf_ =zeros(1e3,2*T+1);
    %
    SDAltGDFed_ = zeros(1e3,2*T+1);
    timeAltGDFed_ =zeros(1e3,2*T+1);
    %
    %SDAltGDFed1_ = zeros(1e3,2*T+1);
    %timeAltGDFed1_ =zeros(1e3,2*T+1);
    %
    SDAltGDMineta1_ = zeros(1e3,2*T+1);
    timeAltGDMineta1_ = zeros(1e3,2*T+1);   
    % --
    strt = 1;
    for i = 1 : numData
        %load("n_"+num2str(n)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_p_"+num2str(p),...
        %      +"_MC_"+num2str(mc) + "_randID_" + num2str(ID) + ".mat")
        %load([num2str(i) +'.mat'])
        load(saveName)
        % --
        SDAltMinParfor_(strt:strt+mc-1,:) = SDAltMinParfor;
        timeAltMinParfor_(strt:strt+mc-1,:) = timeAltMinParfor;
        % 
        %SDAltMinPrvt_(strt:strt+mc-1,:) = SDAltMinPrvt;
        %timeAltMinPrvt_(strt:strt+mc-1,:) = timeAltMinPrvt;
        %
        %SDAltGDFedHalf_(strt:strt+mc-1,:) = SDAltGDFedHalf;
        %timeAltGDFedHalf_(strt:strt+mc-1,:) = timeAltGDFedHalf;        
        %
        SDAltGDFed_(strt:strt+mc-1,:) = SDAltGDFed;
        timeAltGDFed_(strt:strt+mc-1,:) = timeAltGDFed;
        %
        %SDAltGDFed1_(strt:strt+mc-1,:) = SDAltGDFed1;
        %timeAltGDFed1_(strt:strt+mc-1,:) = timeAltGDFed1;
        %
        SDAltGDMineta1_(strt:strt+mc-1,:) = SDAltGDMineta1;
        timeAltGDMineta1_(strt:strt+mc-1,:) = timeAltGDMineta1;
        % -- 
        strt = strt + mc;
    end
        MC = strt; %= legnth of SD/Time arrays
        % -- 
        SDAltMinParfor  = SDAltMinParfor_(1:MC,:);
        timeAltMinParfor = timeAltMinParfor_(1:MC,:);
        % --
        %SDAltMinPrvt = SDAltMinPrvt_(1:MC,:);
        %timeAltMinPrvt = timeAltMinPrvt_(1:MC,:);
        % -- 
        %SDAltGDFedHalf = SDAltGDFedHalf_(1:MC,:);
        %timeAltGDFedHalf = timeAltGDFedHalf_(1:MC,:);
        % -- 
        SDAltGDFed = SDAltGDFed_(1:MC,:);
        timeAltGDFed = timeAltGDFed_(1:MC,:);
        % -- 
        %SDAltGDFed1 = SDAltGDFed1_(1:MC,:);
        %timeAltGDFed1 = timeAltGDFed1_(1:MC,:);
        % -- 
        SDAltGDMineta1 = SDAltGDMineta1_(1:MC,:);
        timeAltGDMineta1 = timeAltGDMineta1_(1:MC,:);
        %---
        % 1
        SDAltMinParfor = sum(SDAltMinParfor,1)/MC;
        timeAltMinParfor = sum(timeAltMinParfor,1)/MC;
        % 3
        %SDAltMinPrvt = sum(SDAltMinPrvt,1)/MC;
        %timeAltMinPrvt = sum(timeAltMinPrvt,1)/MC;
        % 
        %SDAltGDFedHalf = sum(SDAltGDFedHalf,1)/MC;
        %timeAltGDFedHalf = sum(timeAltGDFedHalf,1)/MC;
        % 5
        SDAltGDFed = sum(SDAltGDFed,1)/MC;
        timeAltGDFed = sum(timeAltGDFed,1)/MC;
        % 
        %SDAltGDFed1 = sum(SDAltGDFed1,1)/MC;
        %timeAltGDFed1 = sum(timeAltGDFed1,1)/MC;
        % 7
        SDAltGDMineta1 = sum(SDAltGDMineta1,1)/MC;
        timeAltGDMineta1 = sum(timeAltGDMineta1,1)/MC;
        %------------ Downsample 
        % 1
        SDAltMinParfor = downsample(SDAltMinParfor,1);
        timeAltMinParfor = downsample(timeAltMinParfor,1);
        % 3
        %SDAltMinPrvt = downsample(SDAltMinPrvt,1);
        %timeAltMinPrvt = downsample(timeAltMinPrvt,1);
        % 
        %SDAltGDFedHalf = downsample(SDAltGDFedHalf,4);
        %timeAltGDFedHalf = downsample(timeAltGDFedHalf,4);
        % 5
        SDAltGDFed = downsample(SDAltGDFed,4);
        timeAltGDFed = downsample(timeAltGDFed,4);
        % 
        %SDAltGDFed1 = downsample(SDAltGDFed1,4);
        %timeAltGDFed1 = downsample(timeAltGDFed1,4);
        % 7
        SDAltGDMineta1 = downsample(SDAltGDMineta1,5);
        timeAltGDMineta1 = downsample(timeAltGDMineta1,5);
        %timeSD = 36;
        figure
        semilogy(timeAltGDMineta1(timeAltGDMineta1 <= timeSD ),...
                 SDAltGDMineta1(timeAltGDMineta1 <= timeSD ),'DisplayName', ...
                'AltGDMin(Fed./Prvt.)', ...
                'LineWidth',1.35,'Marker','square','MarkerSize',9, 'Color',"#0000FF")
        hold on
        semilogy(timeAltMinParfor(timeAltMinParfor <= timeSD ),...
                 SDAltMinParfor(timeAltMinParfor <= timeSD),'DisplayName', ...
                'AltMin(Fed./NotPrvt.)', ...
                'LineWidth',1.35,'Marker','x','MarkerSize',9,'Color',"#D95319")
        %semilogy(timeAltMinPrvt(timeAltMinPrvt <= timeSD ),...
        %        SDAltMinPrvt(timeAltMinPrvt <= timeSD),'DisplayName', ...
        %       'AltMin(Fed./Prvt.)', ...
        %       'LineWidth',1.35,'Marker','*','MarkerSize',9,'Color',"#EDB120")
        %-
        %semilogy(timeAltGDFedHalf(timeAltGDFedHalf <= timeSD ),...
        %       SDAltGDFedHalf(timeAltGDFedHalf <= timeSD),'DisplayName', ...
        %      'FactGD (Fed./Prvt. $c = 0.5$)', ...
        %      'LineWidth',1.35,'Marker','diamond','MarkerSize',9,'Color',"0.72,0.27,1.00")
        %-
        semilogy(timeAltGDFed(timeAltGDFed <= timeSD ),...
               SDAltGDFed(timeAltGDFed <= timeSD),'DisplayName', ...
              'FactGD (Fed./Prvt. $c = 0.75$)', ...
              'LineWidth',1.35,'Marker','diamond','MarkerSize',9,'Color',"#7E2F8E")
        %-
        %semilogy(timeAltGDFed1(timeAltGDFed1 <= timeSD ),...
        %       SDAltGDFed1(timeAltGDFed1 <= timeSD),'DisplayName', ...
        %      'FactGD (Fed./Prvt. $c = 1.0$)', ...
        %      'LineWidth',1.35,'Marker','diamond','MarkerSize',9,'Color',"#7E2FEE")
        grid on
        legend('Interpreter','Latex','Fontsize',9.25, 'Location','Northeast')
        yHndl = ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$','Interpreter','Latex','Fontsize',14.0);
        %yHndl.Position(1) = yHndl.Position(1) + 1.25;
        %yHndl.Position(2) = yHndl.Position(2) - 0.05*yHndl.Position(2);
        xlabel('t/seconds','Fontsize',13,'Interpreter','Latex')
        title("n = " + n + ", q = " + q +...
              ", r = " + r + ", p = " + p + ", N $\sim \mathcal{N}$ (0, " + noiseVar_ + "), Wrkrs = " + num2str(numWrkrs),...
               'Interpreter', 'Latex', 'FontSize',14)
        cores = feature('numCores');
        stringTitle = ['Fed_Wrkrs_',num2str(numWrkrs),'Max',num2str(cores),'_MC_',num2str(MC),...
                          '_n_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                           'T_',num2str(T),'_id',num2str(randi(1e3,1))];
        cd(dir)
        exportgraphics(gcf,[stringTitle,'.pdf'])
end
%---
function [SDVals, timeArr] = factGDEta(Xzeros,r,Ustr,T,p,rowIdx,colIdx,Xcol,numWrkrs,Tsvd,eta)
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
    step_const = eta;
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