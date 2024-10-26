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
r = 2;
n = 500;
q = 1000;
% sub sample X with probaility p
p = 0.05;
numWrkrs = 5;
space = 25;
T = 50 + space;
delete(gcp('nocreate'))
parpool(numWrkrs,'IdleTimeout',Inf,'SpmdEnabled',false)
pool = gcp;
MC = 3;
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
%---------------------------------------
pltFedOnly(0, 0, ...
            0, 0, 0,0, ...
            0,0, 0, 0,...
            timeAltGDMineta1,timeAltGDMineta1QR, SDAltGDMineta1, ...
            n, q, r, p, numWrkrs, MC, T,Tsvd)

%saveFigFed(1,dir,21,saveName)
%---
function pltFedOnly(timeAltMinParfor, SDAltMinParfor, ...
                    timeAltMinPrvt, SDAltMinPrvt, timeAltGDFedHalf, SDAltGDFedHalf, ...
                    timeAltGDFed,SDAltGDFed, timeAltGDFed1,SDAltGDFed1,...
                    timeAltGDMineta1, timeAltGDMineta1QR, SDAltGDMineta1,...
                    n, q, r, p, numWrkrs, MC, T,Tsvd)

    %idx1 = find(SDAltMinParfor < 10^-14,1);
    %tidx1 = timeAltMinParfor(idx1);
    %
    %idx3 = find(SDAltMinPrvt < 10^-14,1);
    %tidx3 = timeAltMinPrvt(idx3);
    %
    %idx4 = find(SDAltGDFedHalf  < 10^-15,1);
    %tidx4 = timeAltGDFedHalf(idx4);
    %
    %idx5 = find(SDAltGDFed  < 10^-15,1);
    %tidx5 = timeAltGDFed(idx5);
    %
    %idx6 = find(SDAltGDFed1  < 10^-15,1);
    %tidx6 = timeAltGDFed1(idx6);
    %
    idx7 = find(SDAltGDMineta1 < 10^-14,1);
    tidx7 =  timeAltGDMineta1(idx7);
    %
    timeSD = max([tidx7,tidx7,tidx7,tidx7,tidx7,tidx7]);
    %
    figure;
    %semilogy(timeAltMinParfor(timeAltMinParfor <= timeSD), ...
    %    SDAltMinParfor(timeAltMinParfor <= timeSD), ...
    %    'DisplayName', 'AltMin(Fed./NotPrvt.)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 5);
    %hold on;
    %semilogy(timeAltMinPrvt(timeAltMinPrvt <= timeSD), ...
    %    SDAltMinPrvt(timeAltMinPrvt <= timeSD), ...
    %    'DisplayName', 'AltMin(Fed./Prvt.)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 5);
    %
    %semilogy(timeAltGDFedHalf(timeAltGDFedHalf <= timeSD), ...
    %    SDAltGDFedHalf(timeAltGDFedHalf <= timeSD), ...
    %    'DisplayName', 'FactGD (Fed. $c = 0.5$)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    %
    %semilogy(timeAltGDFed(timeAltGDFed <= timeSD), ...
    %    SDAltGDFed(timeAltGDFed <= timeSD), ...
    %    'DisplayName', 'FactGD (Fed. $c = 0.75$)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    %
    %semilogy(timeAltGDFed1(timeAltGDFed1 <= timeSD), ...
    %    SDAltGDFed1(timeAltGDFed1 <= timeSD), ...
    %    'DisplayName', 'FactGD (Fed. $c = 1.00$)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    % 
    semilogy(timeAltGDMineta1(timeAltGDMineta1 <= timeSD), ...
        SDAltGDMineta1(timeAltGDMineta1 <= timeSD), ...
        'DisplayName', 'AltGDMin(Fed. QR + LS)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 5);
    hold on
    semilogy(timeAltGDMineta1QR(timeAltGDMineta1QR <= timeSD), ...
        SDAltGDMineta1(timeAltGDMineta1QR <= timeSD), ...
        'DisplayName', 'AltGDMin (Fed. QR Only)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 5);
    grid on;
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$', 'Interpreter', 'Latex', 'Fontsize', 15);
    xlabel('Time/seconds', 'FontSize', 11);
    
    title("n = " + n + ", q = " + q +...
        ", r = " + r + ", p = " + p + '.', ...
        'Interpreter', 'Latex', 'FontSize', 14);
    %
    %cores = feature('numCores');
    %stringTitle = ['FedOnly_', num2str(numWrkrs),'_MC_', num2str(MC), ...
    %    '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_p_', num2str(p), ...
    %    'T_', num2str(T),'_Tsvd',num2str(Tsvd),'_id', num2str(randi(1e3, 1))];
    
    %savefig([stringTitle, '.fig']);
end
%---
