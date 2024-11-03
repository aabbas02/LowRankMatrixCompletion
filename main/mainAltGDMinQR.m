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
numWrkrs = 4;
space = 25;
T = 50 + space;
delete(gcp('nocreate'))
parpool(numWrkrs,'IdleTimeout',Inf,'SpmdEnabled',false)
pool = gcp;
MC = 5;
tic 
ID = randi(1e3)+5
ID = 7891359
ID = 2
%--------------------------------------------------------------------------
% generate rank-r X*
U = orth(randn(n,r));
Bstar = randn(r,q);
X  = U*Bstar;
Ustr = U(:,1:r);
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
    % --- AltGDMin
    kAltGDMin = 2;
    eta_c = 1.00;
    [SDAltGDMineta1(mc,:),timeAltGDMineta1(mc,:),timeAltGDMineta1QR(mc,:)] = altGDMin_T(r,eta_c, ...
                                                               Ustr,Xzeros, kAltGDMin*T,p, ...
                                                               rowIdx,Xcol,numWrkrs,Tsvd);
   mc
end
toc
% 7
MC = mc;
SDAltGDMineta1 = sum(SDAltGDMineta1,1)/MC;
timeAltGDMineta1 = sum(timeAltGDMineta1,1)/MC;
timeAltGDMineta1QR = sum(timeAltGDMineta1QR,1)/MC;
%7
SDAltGDMineta1 = downsample(SDAltGDMineta1,5);
timeAltGDMineta1 = downsample(timeAltGDMineta1,5);
timeAltGDMineta1QR = downsample(timeAltGDMineta1QR,5);
figure
timeSD = timeAltGDMineta1(find(SDAltGDMineta1 < 1e-14,1));
semilogy(timeAltGDMineta1(timeAltGDMineta1 <= timeSD ),...
         SDAltGDMineta1(timeAltGDMineta1 <= timeSD ),'DisplayName', ...
        'AltGDMin(Fed./QR + LS)', ...
        'LineWidth',1.35,'Marker','square','MarkerSize',9, 'Color',"#0000FF")
%hold on
grid on
legend('Interpreter','Latex','Fontsize',9.25, 'Location','Northeast')
yHndl = ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$','Interpreter','Latex','Fontsize',14.0);
yHndl.Position(1) = yHndl.Position(1) + 1.25;
yHndl.Position(2) = yHndl.Position(2) - 0.05*yHndl.Position(2);
xlabel('t/seconds','Fontsize',13,'Interpreter','Latex')
title("n = " + n + ", q = " + q +...
      ", r = " + r + ", p = " + p + ", Workers =" + num2str(numWrkrs),...
       'Interpreter', 'Latex', 'FontSize',14)
cores = feature('numCores');
stringTitle = ['Fed_Wrkrs_AltGDMin_Full',num2str(numWrkrs),'Max',num2str(cores),'_MC_',num2str(MC),...
                  '_n_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                   'T_',num2str(T),'_id',num2str(randi(1e3,1))];
exportgraphics(gcf,[stringTitle,'.pdf'])
figure
semilogy(timeAltGDMineta1QR(timeAltGDMineta1QR <= timeSD ),...
         SDAltGDMineta1(timeAltGDMineta1QR <= timeSD ),'DisplayName', ...
        'AltGDMin(Fed./QR Only) ', ...
        'LineWidth',1.35,'Marker','square','MarkerSize',9, 'Color',"#AA00FF")
grid on
legend('Interpreter','Latex','Fontsize',9.25, 'Location','Northeast')
yHndl = ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$','Interpreter','Latex','Fontsize',14.0);
yHndl.Position(1) = yHndl.Position(1) + 1.25;
yHndl.Position(2) = yHndl.Position(2) - 0.05*yHndl.Position(2);
xlabel('t/seconds','Fontsize',13,'Interpreter','Latex')
title("n = " + n + ", q = " + q +...
      ", r = " + r + ", p = " + p + ", Workers =" + num2str(numWrkrs),...
       'Interpreter', 'Latex', 'FontSize',14)
cores = feature('numCores');
stringTitle = ['Fed_Wrkrs_AltGDMin_Only_QR_',num2str(numWrkrs),'Max',num2str(cores),'_MC_',num2str(MC),...
                  '_n_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                   'T_',num2str(T),'_id',num2str(randi(1e3,1))];
cd(dir)
exportgraphics(gcf,[stringTitle,'.pdf'])
