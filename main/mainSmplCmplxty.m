clc
close all
clear all
r = 10;
n = 5000;
q = 5000;
% sub sample X with probaility p
p = 0.015;
numWrkrs = 10;
space = 25;
T = 325 + space;
delete(gcp('nocreate'))
parpool(numWrkrs,'IdleTimeout',Inf,'SpmdEnabled',false)
pool = gcp;
MC = 125;
tic 
ID = randi(1e3);
%--------------------------------------------------------------------------
% generate rank-r X*
U = orth(randn(n,r));
Bstar = randn(r,q);
X  = U*Bstar;
Ustr = U(:,1:r);
for mc = 1 : MC
    % U = orth(randn(n,r));
    % Bstar = randn(r,q);
    % X  = U*Bstar;
    % Ustr = U(:,1:r);
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
    Tsvd = 15;
    %-----------------------------------------------------------------------------------
    % --- AltGD (Federated)
    kAltGD = 2;
    [SDAltGDFed(mc,:), timeAltGDFed(mc,:)] = factGDNew(Xzeros,r,Ustr,T,p,rowIdx,colIdx,Xcol,numWrkrs,Tsvd);
    % --- AltMin (Parfor)
    [SDAltMinParfor(mc,:), timeAltMinParfor(mc,:)] = altMinParfor(Xzeros,r,p, ...
                                          Ustr,T, ...
                                          rowIdx,Xcol,colIdx,Xrow,numWrkrs,Tsvd);                                      
    % --- AltGDMin
    kAltGDMin = 2;
    eta_c = 1.0;
    [SDAltGDMineta1(mc,:),timeAltGDMineta1(mc,:)] = altGDMin(r,eta_c,...
                                                             Ustr,Xzeros, T,p, ...
                                                             rowIdx,Xcol,numWrkrs,Tsvd);
    % if  (mod(mc,5) == 0)
    %     save("n_"+num2str(n)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_p_"+num2str(p)+"_MC_"+num2str(mc) + "_randID_" + num2str(ID) + ".mat",...
    %                            "SDAltGDMineta1","timeAltGDMineta1",...
    %                            "SDAltMinParfor","timeAltMinParfor",...
    %                            "SDAltGDFed","timeAltGDFed",...
    %                            "n","p","q","r","numWrkrs","mc","T",...
    %                            "numWrkrs");             %"SDAltGDMin","timeAltGDMin",...
                   
    %end
    mc
end
toc
% 1
MC = mc;
SDAltMinParfor = sum(SDAltMinParfor,1)/MC;
timeAltMinParfor = sum(timeAltMinParfor,1)/MC;
% 4
SDAltGDFed = sum(SDAltGDFed,1)/MC;
timeAltGDFed = sum(timeAltGDFed,1)/MC;
%
SDAltGDMineta1 = sum(SDAltGDMineta1,1)/MC;
timeAltGDMineta1 = sum(timeAltGDMineta1,1)/MC;
%--- Plot Subspace Distance against Time Figure
idx1 = find(SDAltMinParfor < 10^-14,1);
tidx1 = timeAltMinParfor(idx1);

idx5 = find(SDAltGDFed  < 10^-15,1);
tidx5 = timeAltGDFed(idx5);


idx7 = find(SDAltGDMineta1 < 10^-14,1);
tidx7 =  timeAltGDMineta1(idx7);


timeSD = max([tidx1,tidx5,tidx7]);
if (isempty(timeSD))
    timeSD = 100;
end

figure
semilogy(timeAltMinParfor(timeAltMinParfor <= timeSD ),...
         SDAltMinParfor(timeAltMinParfor <= timeSD),'DisplayName', ...
        'AltMin(Fed./NotPrvt.)', ...
        'LineWidth',1.45,'Marker','o','MarkerSize',5)
hold on
semilogy(timeAltGDFed(timeAltGDFed <= timeSD ),...
        SDAltGDFed(timeAltGDFed <= timeSD),'DisplayName', ...
       'Fact GD(Fed. Prvt.)', ...
       'LineWidth',1.45,'Marker','x','MarkerSize',5)
semilogy(timeAltGDMineta1(timeAltGDMineta1 <= timeSD ),...
         SDAltGDMineta1(timeAltGDMineta1 <= timeSD ),'DisplayName', ...
        'AltGDMin(Fed.)', ...
        'LineWidth',1.45,'Marker','square','MarkerSize',5)
grid on

legend('Interpreter','Latex','Fontsize',9)
ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$','Interpreter','Latex','Fontsize',15)
xlabel('Time/seconds',FontSize=11)
title("n = " + n + ", q = " + q +...
      ", r = " + r + ", p = " + p + '.',...
       'Interpreter', 'Latex', 'FontSize',14)
cores = feature('numCores');
stringTitle = ['notSparse_Wrkrs_',num2str(numWrkrs),'Max',num2str(cores),'_MC_',num2str(MC),...
                  '_n_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                   'T_',num2str(T),'_id',num2str(randi(1e3,1))];

savefig([stringTitle,'.fig']);
cd(dir)    
savefig([stringTitle, '.fig']);
%exportgraphics(gcf,[stringTitle,'.pdf'])


