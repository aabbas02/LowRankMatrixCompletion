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
r = 10;
n = 5000;
q = 5000;
% sub sample X with probaility p
p = 0.1;
numWrkrs = 10;
space = 25;
T = 25 + space;
delete(gcp('nocreate'))
parpool(numWrkrs,'IdleTimeout',Inf,'SpmdEnabled',false)
pool = gcp;
MC = 1;
tic 
% generate rank-r X*
Ustr = orth(randn(n,r));
Bstar = randn(r,q);
X  = Ustr*Bstar;
for mc = 1 : MC
    %Instantiate Y = X_Omega
    % Make cell arrays for row and column indices
    %-----------------------------------------------------------------------------------
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
