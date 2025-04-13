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
q = 5000;
% sub sample X with probaility p
p = 0.05;
numWrkrs = 10;
space = 25;
T = 25 + space;
delete(gcp('nocreate'))
parpool(numWrkrs,'IdleTimeout',Inf,'SpmdEnabled',false)
pool = gcp;
MC = 1;
tic 
ID = randi(1e3)+5
ID = 7891359
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
end
parfor j = 1 : n
    colIdx{j} = col(row==j)
    Xrow{j} = X(j,colIdx{j})';    
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
    end
    parfor j = 1 : n
       colIdx{j} = col(row==j)
       Xrow{j} = X(j,colIdx{j})';        
    end        
    %-----------------------------------------------------------------------------------
    % --- AltGD (Federated)
    kAltGD = 2;
    [SDAltGDFed(mc,:), timeAltGDFed(mc,:)] = factGDNew(Xzeros,r,Ustr,kAltGD*T,p,...
                                                       rowIdx,colIdx,Xcol,numWrkrs,Tsvd);
    % --- AltMin (Parfor)
    [SDAltMinParfor(mc,:), timeAltMinParfor(mc,:)] = altMinParfor_T(Xzeros,r,p, ...
									                              Ustr,T, ...
									                              rowIdx,Xcol,colIdx,Xrow,numWrkrs,Tsvd);
    % --- AltMin (Cntrl.)
    [SDAltMinCntrl(mc,:), timeAltMinCntrl(mc,:)] = altMinCntrl_T(Xzeros,r,p, ...
                                                                 Ustr,T, ...
                                                                 rowIdx,Xcol,colIdx,Xrow,numWrkrs,Tsvd);
    % --- AltMin (Fed./Prvt.)
    T_inner =  10;
    [SDAltMinPrvt(mc,:), timeAltMinPrvt(mc,:)]  = altMinPrvt_T(Xzeros,r,p, ...
                                                               Ustr,T, ...
                                                               rowIdx,Xcol,T_inner,numWrkrs,Tsvd);

    % --- AltGD (Cntrl.)
    kAltGD = 2;
    [SDAltGD(mc,:), timeAltGD(mc,:)] = AltGD(Xzeros,r,Ustr,kAltGD*T,p,idxC,numWrkrs,Tsvd);                                                     
    % --- AltGDMin
    kAltGDMin = 2;
    eta_c = 1.00;
    [SDAltGDMineta1(mc,:),timeAltGDMineta1(mc,:)] = altGDMin_T(r,eta_c, ...
                                                               Ustr,Xzeros, kAltGDMin*T,p, ...
                                                               rowIdx,Xcol,numWrkrs,Tsvd);
    % --- AltGDMin(Cntrl.)
    [rowSort, colSort, ~] = find(Xzeros);
    idx = sub2ind([n,q],rowSort,colSort);
    Xvec = X(idx);
    [SDAltGDMinCntrl(mc,:),timeAltGDMinCntrl(mc,:)] = altGDMinCntrl_T(Xzeros,Xvec,r,p,...
                                                                      rowSort,colSort,rowIdx,Xcol,...
                                                                      T,Ustr,numWrkrs,Tsvd);

    if  (mod(mc,5) == 0)
        save("n_"+num2str(n)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_p_"+num2str(p)+"_MC_"+num2str(mc) + "_randID_" + num2str(ID) + ".mat",...
                               "SDAltGDMineta1","timeAltGDMineta1",...
                               "SDAltMinParfor","timeAltMinParfor",...
                               "SDAltMinCntrl","timeAltMinCntrl",...
                               "SDAltGDMinCntrl","timeAltGDMinCntrl",...
                               "SDAltMinPrvt","timeAltMinPrvt",...
                               "SDAltGDFed","timeAltGDFed",...
                               "SDAltGD","timeAltGD",...
                               "n","p","q","r","numWrkrs","mc","T",...
                               "numWrkrs");             %"SDAltGDMin","timeAltGDMin",...
                   
   end
    mc
end
toc
% 1
MC = mc;
SDAltMinParfor = sum(SDAltMinParfor,1)/MC;
timeAltMinParfor = sum(timeAltMinParfor,1)/MC;
% 2
SDAltMinCntrl = sum(SDAltMinCntrl,1)/MC;
timeAltMinCntrl = sum(timeAltMinCntrl,1)/MC;
% 3
SDAltMinPrvt = sum(SDAltMinPrvt,1)/MC;
timeAltMinPrvt = sum(timeAltMinPrvt,1)/MC;
% 4
SDAltGDFed = sum(SDAltGDFed,1)/MC;
timeAltGDFed = sum(timeAltGDFed,1)/MC;
% 5
SDAltGD = sum(SDAltGD,1)/MC;
timeAltGD = sum(timeAltGD,1)/MC;
% 6
%SDAltGDMin = sum(SDAltGDMin,1)/MC;
%timeAltGDMin = sum(timeAltGDMin,1)/MC;
% 7
SDAltGDMineta1 = sum(SDAltGDMineta1,1)/MC;
timeAltGDMineta1 = sum(timeAltGDMineta1,1)/MC;
% 8
SDAltGDMinCntrl = sum(SDAltGDMinCntrl,1)/MC;
timeAltGDMinCntrl = sum(timeAltGDMinCntrl,1)/MC;
% --- Plot Subspace Distance against Iteration Figure
%plotErrvsIter(space, T, n, q, r, p,...
%              SDAltMinCntrl, SDAltMinParfor, SDAltMinPrvt, SDAltGDFed,...
%              SDAltGD, SDAltGDMineta1)
%--- Plot Subspace Distance against Time Figure
plotAndSaveFigureFull(timeAltMinCntrl, SDAltMinCntrl, timeAltMinParfor, SDAltMinParfor, ...
    timeAltMinPrvt, SDAltMinPrvt, timeAltGD, SDAltGD, timeAltGDFed, SDAltGDFed, ...
    timeAltGDMineta1, SDAltGDMineta1, timeAltGDMinCntrl, SDAltGDMinCntrl, ...
    n, q, r, p, numWrkrs, MC, T,Tsvd)

stringTitle = ['Lab_Wrkrs', num2str(numWrkrs),'_MC_', num2str(MC), ...
               '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_p_', num2str(p), ...
               'T_', num2str(T),'_Tsvd',num2str(Tsvd),'_id', num2str(randi(1e3, 1))];

savefig([stringTitle, '.fig']);
%exportgraphics(gcf,[stringTitle,'.pdf'])
