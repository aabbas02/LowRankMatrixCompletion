clc
close all
clear all
dir = pwd;
% % For linux, replace '\' with '/'
% cd ('C:\Users\aabbasi1\OneDrive - Iowa State University\Desktop\FLRMC\altGDMin2ndPaper')
% addpath(genpath('.\functions'));
% numData = 4;
% saveFig(numData,dir)
% clc
% close all
% clear all
% dir = pwd;
% % For linux, replace '\' with '/'
% cd ('C:\Users\aabbasi1\OneDrive - Iowa State University\Desktop\FLRMC\altGDMin2ndPaper')
% addpath(genpath('.\functions'));
% numData = 4;
% saveFigCntrl(numData,dir)
% clc
% close all
% clear all
% dir = pwd;
% For linux, replace '\' with '/'
%cd ('C:\Users\aabbasi1\OneDrive - Iowa State University\Desktop\FLRMC\altGDMin2ndPaper')
numData = 3;
saveFigFed(numData,dir,15)
close all
clc


function saveFigFed(numData,dir,timeSD)
    cd (dir)
    clc
    close all
    load("1.mat");
    %T = 0; % You need to define the value of T, as it is used in the code
    numWrkrs = 10; % Define the number of workers
    % --
    SDAltMinParfor_ = zeros(1e3,T+1);
    timeAltMinParfor_ = zeros(1e3,T+1);
    SDAltMinPrvt_ = zeros(1e3,T+1);
    timeAltMinPrvt_ = zeros(1e3,T+1);
    SDAltGDFed_ = zeros(1e3,2*T+1);
    timeAltGDFed_ =zeros(1e3,2*T+1);
    SDAltGDMineta1_ = zeros(1e3,2*T+1);
    timeAltGDMineta1_ = zeros(1e3,2*T+1);   
    % --
    strt = 1;
    for i = 1 : numData
        load([num2str(i+4) +'.mat'])
        % --
        SDAltMinParfor_(strt:strt+mc-1,:) = SDAltMinParfor;
        timeAltMinParfor_(strt:strt+mc-1,:) = timeAltMinParfor;
        SDAltMinPrvt_(strt:strt+mc-1,:) = SDAltMinPrvt;
        timeAltMinPrvt_(strt:strt+mc-1,:) = timeAltMinPrvt;
        SDAltGDFedHalf_(strt:strt+mc-1,:) = SDAltGDFedHalf;
        timeAltGDFedHalf_(strt:strt+mc-1,:) = timeAltGDFedHalf;
        SDAltGDFed_(strt:strt+mc-1,:) = SDAltGDFed;
        timeAltGDFed_(strt:strt+mc-1,:) = timeAltGDFed;
        SDAltGDFed1_(strt:strt+mc-1,:) = SDAltGDFed1;
        timeAltGDFed1_(strt:strt+mc-1,:) = timeAltGDFed1;
        SDAltGDMineta1_(strt:strt+mc-1,:) = SDAltGDMineta1;
        timeAltGDMineta1_(strt:strt+mc-1,:) = timeAltGDMineta1;
        % -- 
        strt = strt + mc-1;
    end
        MC = strt; %= legnth of SD/Time arrays
        % -- 
        SDAltMinParfor  = SDAltMinParfor_(1:MC,:);
        timeAltMinParfor = timeAltMinParfor_(1:MC,:);
        % --
        SDAltMinPrvt = SDAltMinPrvt_(1:MC,:);
        timeAltMinPrvt = timeAltMinPrvt_(1:MC,:);
        % -- 
        SDAltGDFedHalf = SDAltGDFedHalf_(1:MC,:);
        timeAltGDFedHalf = timeAltGDFedHalf_(1:MC,:);
        % -- 
        SDAltGDFed = SDAltGDFed_(1:MC,:);
        timeAltGDFed = timeAltGDFed_(1:MC,:);
        % -- 
        SDAltGDFed1 = SDAltGDFed1_(1:MC,:);
        timeAltGDFed1 = timeAltGDFed1_(1:MC,:);
        % -- 
        SDAltGDMineta1 = SDAltGDMineta1_(1:MC,:);
        timeAltGDMineta1 = timeAltGDMineta1_(1:MC,:);
        %---------------------------------------------
        % 1
        SDAltMinParfor = sum(SDAltMinParfor,1)/MC;
        timeAltMinParfor = sum(timeAltMinParfor,1)/MC;
        % 3
        SDAltMinPrvt = sum(SDAltMinPrvt,1)/MC;
        timeAltMinPrvt = sum(timeAltMinPrvt,1)/MC;
        % -- 
        SDAltGDFedHalf = sum(SDAltGDFedHalf,1)/MC;
        timeAltGDFedHalf = sum(timeAltGDFedHalf,1)/MC;        
        % 5
        SDAltGDFed = sum(SDAltGDFed,1)/MC;
        timeAltGDFed = sum(timeAltGDFed,1)/MC;
        % -- 
        SDAltGDFed1 = sum(SDAltGDFed1,1)/MC;
        timeAltGDFed1 = sum(timeAltGDFed1,1)/MC;
        % 7
        SDAltGDMineta1 = sum(SDAltGDMineta1,1)/MC;
        timeAltGDMineta1 = sum(timeAltGDMineta1,1)/MC;
        %------------ Downsample 
        % 1
        SDAltMinParfor = downsample(SDAltMinParfor,1);
        timeAltMinParfor = downsample(timeAltMinParfor,1);
        % 2
        % 3
        SDAltMinPrvt = downsample(SDAltMinPrvt,1);
        timeAltMinPrvt = downsample(timeAltMinPrvt,1);
        % --
        SDAltGDFedHalf = downsample(SDAltGDFedHalf,4);
        timeAltGDFedHalf = downsample(timeAltGDFedHalf,4);        
        % 5
        SDAltGDFed = downsample(SDAltGDFed,4);
        timeAltGDFed = downsample(timeAltGDFed,4);
        % -- 
        SDAltGDFed1 = downsample(SDAltGDFed1,4);
        timeAltGDFed1 = downsample(timeAltGDFed1,4);
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
        semilogy(timeAltMinPrvt(timeAltMinPrvt <= timeSD ),...
                SDAltMinPrvt(timeAltMinPrvt <= timeSD),'DisplayName', ...
               'AltMin(Fed./Prvt.)', ...
               'LineWidth',1.35,'Marker','*','MarkerSize',9,'Color',"#EDB120")
        %-
        semilogy(timeAltGDFedHalf(timeAltGDFedHalf <= timeSD ),...
               SDAltGDFedHalf(timeAltGDFedHalf <= timeSD),'DisplayName', ...
              'FactGD(Fed./Prvt.$c = 0.5$)', ...
              'LineWidth',1.15,'Marker','diamond','MarkerSize',9,'Color',[0.72,0.27,1.00,0.5])
        %-
        semilogy(timeAltGDFed(timeAltGDFed <= timeSD ),...
               SDAltGDFed(timeAltGDFed <= timeSD),'DisplayName', ...
              'FactGD(Fed./Prvt.$c = 3/4$)', ...
              'LineWidth',1.35,'Marker','diamond','MarkerSize',9,'Color',"#7E2F8E")
        %-
        semilogy(timeAltGDFed1(timeAltGDFed1 <= timeSD ),...
               SDAltGDFed1(timeAltGDFed1 <= timeSD),'DisplayName', ...
              'FactGD(Fed./Prvt.$c = 1.0$)', ...
              'LineWidth',1.15,'Marker','diamond','MarkerSize',9,'Color',"#7E2FEE")        
        grid on
        legend('Interpreter','Latex','Fontsize',8.25','Location','Southwest')%
        yHndl = ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$','Interpreter','Latex','Fontsize',14.0);
        yHndl.Position(1) = yHndl.Position(1) + 1.25;
        yHndl.Position(2) = yHndl.Position(2) - 0.05*yHndl.Position(2);
        xlabel('t/seconds','Fontsize',13,'Interpreter','Latex')
        title("n = " + n + ", q = " + q +...
              ", r = " + r + ", p = " + p + ", Workers =" + 10,...
               'Interpreter', 'Latex', 'FontSize',14)
        cores = feature('numCores');
        stringTitle = ['Fed_Wrkrs_',num2str(numWrkrs),'Max',num2str(cores),'_MC_',num2str(MC),...
                          '_n_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                           'T_',num2str(T),'_id',num2str(randi(1e3,1))];
        cd(dir)
        exportgraphics(gcf,[stringTitle,'.pdf'])
end
