clc
close all
clear all
load("1.mat");
SDAltMinParfor_ = zeros(1e3,T+1);
timeAltMinParfor_ = zeros(1e3,T+1);

SDAltMinCntrl_ = zeros(1e3,T+1);
timeAltMinCntrl_ = zeros(1e3,T+1);

SDAltMinPrvt_ = zeros(1e3,T+1);
timeAltMinPrvt_ = zeros(1e3,T+1);

SDAltGD_ = zeros(1e3,2*T+1);
timeAltGD_ =zeros(1e3,2*T+1);

SDAltGDFed_ = zeros(1e3,2*T+1);
timeAltGDFed_ =zeros(1e3,2*T+1);


%SDProjGD_ = zeros(1e3,2*T+1);
%timeProjGD_ = zeros(1e3,2*T+1);

%SDAltGDMin_ =zeros(1e3,2*T+1);
%timeAltGDMin_ = zeros(1e3,2*T+1);

SDAltGDMineta1_ = zeros(1e3,2*T+1);
timeAltGDMineta1_ = zeros(1e3,2*T+1);

SDAltGDMinCntrl_ = zeros(1e3,T+1);
timeAltGDMinCntrl_ = zeros(1e3,T+1);

strt = 1;
for i = 1 : 4
    load([num2str(i) +'.mat'])

    SDAltMinParfor_(strt:strt+mc-1,:) = SDAltMinParfor;
    timeAltMinParfor_(strt:strt+mc-1,:) = timeAltMinParfor;

    SDAltMinCntrl_(strt:strt+mc-1,:) = SDAltMinCntrl;
    timeAltMinCntrl_(strt:strt+mc-1,:) = timeAltMinCntrl;

    SDAltMinPrvt_(strt:strt+mc-1,:) = SDAltMinPrvt;
    timeAltMinPrvt_(strt:strt+mc-1,:) = timeAltMinPrvt;

    %SDProjGD_(strt:strt+mc-1,:) = SDProjGD;
    %timeProjGD_(strt:strt+mc-1,:) = timeProjGD;

    %SDAltGDMin_(strt:strt+mc-1,:) = SDAltGDMin;
    %timeAltGDMin_(strt:strt+mc-1,:) = timeAltGDMin;

    SDAltGD_(strt:strt+mc-1,:) = SDAltGD;
    timeAltGD_(strt:strt+mc-1,:) = timeAltGD;

    SDAltGDFed_(strt:strt+mc-1,:) = SDAltGDFed;
    timeAltGDFed_(strt:strt+mc-1,:) = timeAltGDFed;


    SDAltGDMineta1_(strt:strt+mc-1,:) = SDAltGDMineta1;
    timeAltGDMineta1_(strt:strt+mc-1,:) = timeAltGDMineta1;

    SDAltGDMinCntrl_(strt:strt+mc-1,:) = SDAltGDMinCntrl;
    timeAltGDMinCntrl_(strt:strt+mc-1,:) = timeAltGDMinCntrl;

    strt = strt + mc;
end

MC = strt; %= legnth of SD/Time arrays

SDAltMinParfor  = SDAltMinParfor_(1:MC,:);
timeAltMinParfor = timeAltMinParfor_(1:MC,:);

SDAltMinCntrl = SDAltMinCntrl_(1:MC,:);
timeAltMinCntrl = timeAltMinCntrl_(1:MC,:);

SDAltMinPrvt = SDAltMinPrvt_(1:MC,:);
timeAltMinPrvt = timeAltMinPrvt_(1:MC,:);

SDAltGD = SDAltGD_(1:MC,:);
timeAltGD = timeAltGD_(1:MC,:);

SDAltGDFed = SDAltGDFed_(1:MC,:);
timeAltGDFed = timeAltGDFed_(1:MC,:);


%SDProjGD = SDProjGD_(1:MC,:);
%timeProjGD = timeProjGD_(1:MC,:);

%SDAltGDMin = SDAltGDMin_(1:MC,:);
%timeAltGDMin = timeAltGDMin_(1:MC,:);

SDAltGDMineta1 = SDAltGDMineta1_(1:MC,:);
timeAltGDMineta1 = timeAltGDMineta1_(1:MC,:);

SDAltGDMinCntrl = SDAltGDMinCntrl_(1:MC,:);
timeAltGDMinCntrl = timeAltGDMinCntrl_(1:MC,:);
%---------------------------------------------
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

SDAltGDFed = sum(SDAltGDFed,1)/MC;
timeAltGDFed = sum(timeAltGDFed,1)/MC;

% 5
%SDProjGD = sum(SDProjGD,1)/MC;
%timeProjGD = sum(timeProjGD,1)/MC;
% 6
%SDAltGDMin = sum(SDAltGDMin,1)/MC;
%timeAltGDMin = sum(timeAltGDMin,1)/MC;
% 7
SDAltGDMineta1 = sum(SDAltGDMineta1,1)/MC;
timeAltGDMineta1 = sum(timeAltGDMineta1,1)/MC;
% 8
SDAltGDMinCntrl = sum(SDAltGDMinCntrl,1)/MC;
timeAltGDMinCntrl = sum(timeAltGDMinCntrl,1)/MC;
%----------------------------------------------------------------------

%------------ Downsample 
% 1
SDAltMinParfor = downsample(SDAltMinParfor,1);
timeAltMinParfor = downsample(timeAltMinParfor,1);
% 2
SDAltMinCntrl = downsample(SDAltMinCntrl,1);
timeAltMinCntrl = downsample(timeAltMinCntrl,1);
% 3
SDAltMinPrvt = downsample(SDAltMinPrvt,1);
timeAltMinPrvt = downsample(timeAltMinPrvt,1);
% 4
SDAltGD = downsample(SDAltGD,4);
timeAltGD = downsample(timeAltGD,4);

SDAltGDFed = downsample(SDAltGDFed,4);
timeAltGDFed = downsample(timeAltGDFed,4);

% 5
%SDProjGD = downsample(SDProjGD,1);
%timeProjGD = downsample(timeProjGD,1);
% 6
%SDAltGDMin = downsample(SDAltGDMin,5);
%timeAltGDMin = downsample(timeAltGDMin,5);
% 7
SDAltGDMineta1 = downsample(SDAltGDMineta1,5);
timeAltGDMineta1 = downsample(timeAltGDMineta1,5);
% 8
SDAltGDMinCntrl = downsample(SDAltGDMinCntrl,5);
timeAltGDMinCntrl = downsample(timeAltGDMinCntrl,5);


% --- Plot Subspace Distance against Iteration Figure
%plotErrvsIter(space, T, n, q, r, p,...
%               SDAltMinCntrl, SDAltMinParfor, SDAltMinPrvt, SDAltGD,...
%               SDProjGD, SDAltGDMin, SDAltGDMineta1)
%--- Plot Subspace Distance against Time Figure


idx1 = find(SDAltMinParfor < 10^-14,1);
tidx1 = timeAltMinParfor(idx1);

idx2 = find(SDAltMinCntrl < 10^-14,1);
tidx2 = timeAltMinCntrl(idx2);

idx3 = find(SDAltMinPrvt < 10^-14,1);
tidx3 = timeAltMinPrvt(idx3);

idx4 =  find(SDAltGD < 10^-14,1);
tidx4 = timeAltGD(idx4);

%idx5 = find(SDProjGD < 10^-14,1);
%tidx5 =  timeProjGD(idx5);

%idx6 = find(SDAltGDMin < 10^-14,1);
%tidx6 =  timeAltGDMin(idx6);

idx7 = find(SDAltGDMineta1 < 10^-14,1);
tidx7 =  timeAltGDMineta1(idx7);

idx8 = find(SDAltGDMinCntrl < 10^-14,1);
tidx8 =  timeAltGDMinCntrl(idx8);

%timeSD = max([tidx1,tidx2,tidx3,tidx4,tidx5,tidx6,tidx7,tidx8]);
%timeSD = max([tidx1,tidx2,tidx3,tidx4,tidx5,tidx6,tidx7,tidx8]);
%timeSD = max([tidx1,tidx2,tidx4,tidx6,tidx7,tidx8]);
%timeSD = max([tidx1,tidx7]);
timeSD = 30;


figure

% semilogy(timeAltGDMineta1(timeAltGDMineta1 <= timeSD ),...
%          SDAltGDMineta1(timeAltGDMineta1 <= timeSD ),'DisplayName', ...
%         'AltGDMin(Fed.) $\eta = 1.0$', ...
%         'LineWidth',1.35,'Marker','square','MarkerSize',9, 'Color',"#0000FF")
% 
% 
% hold on
% 
% semilogy(timeAltGDMin(timeAltGDMin <= timeSD ),...
%         SDAltGDMin(timeAltGDMin <= timeSD ),'DisplayName', ...
%        'AltGDMin(Fed.) $\eta = 3/4$', ...
%        'LineWidth',1.35,'Marker','square','MarkerSize',9,'Color',"#750AFF")
% 
% semilogy(timeAltGDMinCntrl(timeAltGDMinCntrl <= timeSD ),...
%         SDAltGDMinCntrl(timeAltGDMinCntrl <= timeSD ),'DisplayName', ...
%        'AltGDMin(Cntrl.)', ...
%        'LineWidth',1.35,'Marker','square','MarkerSize',9,'Color',"#0072BD")
% 
% semilogy(timeAltMinParfor(timeAltMinParfor <= timeSD ),...
%          SDAltMinParfor(timeAltMinParfor <= timeSD),'DisplayName', ...
%         'AltMin(Fed./NotPrvt.)', ...
%         'LineWidth',1.35,'Marker','o','MarkerSize',9,'Color',"#D95319")
% 
% semilogy(timeAltMinPrvt(timeAltMinPrvt <= timeSD ),...
%         SDAltMinPrvt(timeAltMinPrvt <= timeSD),'DisplayName', ...
%        'AltMin(Fed./Prvt.)', ...
%        'LineWidth',1.35,'Marker','o','MarkerSize',9,'Color',"#EDB120")
% 
% semilogy(timeAltMinCntrl(timeAltMinCntrl <= timeSD),...
%        SDAltMinCntrl(timeAltMinCntrl <= timeSD),...
%        'DisplayName', ...
%        'AltMin(Cntrl.)', ...
%        'LineWidth',1.35,'Marker','o','MarkerSize',9,'Color',"#A2142F")
% 
% semilogy(timeAltGD(timeAltGD <= timeSD ),...
%         SDAltGD(timeAltGD <= timeSD),'DisplayName', ...
%        'AltGD', ...
%        'LineWidth',1.35,'Marker','x','MarkerSize',9,'Color',"#77AC30"	)
% 
% semilogy(timeProjGD(timeProjGD <= timeSD ),...
%         SDProjGD(timeProjGD <= timeSD),'DisplayName', ...
%        'ProjGD', ...
%        'LineWidth',1.35,'Marker','+','MarkerSize',9,'Color',"#7E2F8E")

semilogy(timeAltGDMineta1(timeAltGDMineta1 <= timeSD ),...
         SDAltGDMineta1(timeAltGDMineta1 <= timeSD ),'DisplayName', ...
        'AltGDMin(Fed.)', ...
        'LineWidth',1.35,'Marker','square','MarkerSize',9, 'Color',"#0000FF")


hold on

%semilogy(timeAltGDMin(timeAltGDMin <= timeSD ),...
%        SDAltGDMin(timeAltGDMin <= timeSD ),'DisplayName', ...
%       'AltGDMin(Fed.) $\eta = 3/4$', ...
%       'LineWidth',1.35,'Marker','o','MarkerSize',9,'Color',"#750AFF")

semilogy(timeAltGDMinCntrl(timeAltGDMinCntrl <= timeSD ),...
        SDAltGDMinCntrl(timeAltGDMinCntrl <= timeSD ),'DisplayName', ...
       'AltGDMin(Cntrl.)', ...
       'LineWidth',1.35,'Marker','^','MarkerSize',9,'Color',"#0072BD")

semilogy(timeAltMinParfor(timeAltMinParfor <= timeSD ),...
         SDAltMinParfor(timeAltMinParfor <= timeSD),'DisplayName', ...
        'AltMin(Fed./NotPrvt.)', ...
        'LineWidth',1.35,'Marker','x','MarkerSize',9,'Color',"#D95319")

semilogy(timeAltMinPrvt(timeAltMinPrvt <= timeSD ),...
        SDAltMinPrvt(timeAltMinPrvt <= timeSD),'DisplayName', ...
       'AltMin(Fed./Prvt.)', ...
       'LineWidth',1.35,'Marker','*','MarkerSize',9,'Color',"#EDB120")

semilogy(timeAltMinCntrl(timeAltMinCntrl <= timeSD),...
       SDAltMinCntrl(timeAltMinCntrl <= timeSD),...
       'DisplayName', ...
       'AltMin(Cntrl.)', ...
       'LineWidth',1.35,'Marker','hexagram','MarkerSize',9,'Color',"#A2142F")

semilogy(timeAltGD(timeAltGD <= timeSD ),...
        SDAltGD(timeAltGD <= timeSD),'DisplayName', ...
       'AltGD (Cntrl.)', ...
       'LineWidth',1.35,'Marker','square','MarkerSize',9,'Color',"#77AC30"	)

semilogy(timeAltGDFed(timeAltGDFed <= timeSD ),...
       SDAltGDFed(timeAltGDFed <= timeSD),'DisplayName', ...
      'AltGD (Fed.)', ...
      'LineWidth',1.35,'Marker','diamond','MarkerSize',9,'Color',"#7E2F8E")

grid on

legend('Interpreter','Latex','Fontsize',9.25, 'Location','Northeast')
yHndl = ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$','Interpreter','Latex','Fontsize',14.0);
yHndl.Position(1) = yHndl.Position(1) + 1.25;
yHndl.Position(2) = yHndl.Position(2) - 0.05*yHndl.Position(2);
xlabel('t/seconds','Fontsize',13,'Interpreter','Latex')
title("n = " + n + ", q = " + q +...
      ", r = " + r + ", p = " + p + ", Workers =" + 10,...
       'Interpreter', 'Latex', 'FontSize',14)
cores = feature('numCores');
stringTitle = ['Avg_Truncate_',num2str(numWrkrs),'Max',num2str(cores),'_MC_',num2str(MC),...
                  '_n_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                   'T_',num2str(T),'_id',num2str(randi(1e3,1))];

%savefig([stringTitle,'.fig']);
%-------------------------------------------------------------------------
 %timeSD = max([tidx1,tidx7]);
 %ax2 = axes('Position',[0.63 0.65 0.27 0.27]);
 
% semilogy(timeAltGDMineta1(timeAltGDMineta1 <= timeSD ),...
%          SDAltGDMineta1(timeAltGDMineta1 <= timeSD ),'DisplayName', ...
%         'AltGDMin(Fed.) $\eta = 1.0$', ...
%         'LineWidth',1.05,'Marker','square','MarkerSize',7.5,'Color',"#0000FF")
 
% hold on
% grid on
% semilogy(timeAltMinParfor(timeAltMinParfor <= timeSD ),...
 %         SDAltMinParfor(timeAltMinParfor <= timeSD),'DisplayName', ...
 %        'AltMin(Fed./NotPrvt.)', ...
 %        'LineWidth',1.05,'Marker','x','MarkerSize',7.5,'Color',"#D95319")


exportgraphics(gcf,[stringTitle,'.pdf'])


