function saveFigCntrl(numData,dir,numIter)
    cd (dir)
    clc
    close all
    load("1.mat");
    %T = 0; % You need to define the value of T, as it is used in the code
    numWrkrs = 10; % Define the number of workers
    % --
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
    SDAltGDMineta1_ = zeros(1e3,2*T+1);
    timeAltGDMineta1_ = zeros(1e3,2*T+1);   
    SDAltGDMinCntrl_ = zeros(1e3,T+1);
    timeAltGDMinCntrl_ = zeros(1e3,T+1);
    % --
    strt = 1;
    for i = 1 : numData
        load([num2str(i) +'.mat'])
        % --
        SDAltMinParfor_(strt:strt+mc-1,:) = SDAltMinParfor;
        timeAltMinParfor_(strt:strt+mc-1,:) = timeAltMinParfor;
        SDAltMinCntrl_(strt:strt+mc-1,:) = SDAltMinCntrl;
        timeAltMinCntrl_(strt:strt+mc-1,:) = timeAltMinCntrl;
        SDAltMinPrvt_(strt:strt+mc-1,:) = SDAltMinPrvt;
        timeAltMinPrvt_(strt:strt+mc-1,:) = timeAltMinPrvt;
        SDAltGD_(strt:strt+mc-1,:) = SDAltGD;
        timeAltGD_(strt:strt+mc-1,:) = timeAltGD;
        SDAltGDFed_(strt:strt+mc-1,:) = SDAltGDFed;
        timeAltGDFed_(strt:strt+mc-1,:) = timeAltGDFed;
        SDAltGDMineta1_(strt:strt+mc-1,:) = SDAltGDMineta1;
        timeAltGDMineta1_(strt:strt+mc-1,:) = timeAltGDMineta1;
        SDAltGDMinCntrl_(strt:strt+mc-1,:) = SDAltGDMinCntrl;
        timeAltGDMinCntrl_(strt:strt+mc-1,:) = timeAltGDMinCntrl;
        % -- 
        strt = strt + mc;
    end
        MC = strt; %= legnth of SD/Time arrays
        % -- 
        SDAltMinParfor  = SDAltMinParfor_(1:MC,:);
        timeAltMinParfor = timeAltMinParfor_(1:MC,:);
        % -- 
        SDAltMinCntrl = SDAltMinCntrl_(1:MC,:);
        timeAltMinCntrl = timeAltMinCntrl_(1:MC,:);
        % --
        SDAltMinPrvt = SDAltMinPrvt_(1:MC,:);
        timeAltMinPrvt = timeAltMinPrvt_(1:MC,:);
        % --
        SDAltGD = SDAltGD_(1:MC,:);
        timeAltGD = timeAltGD_(1:MC,:);
        % -- 
        SDAltGDFed = SDAltGDFed_(1:MC,:);
        timeAltGDFed = timeAltGDFed_(1:MC,:);
        % -- 
        SDAltGDMineta1 = SDAltGDMineta1_(1:MC,:);
        timeAltGDMineta1 = timeAltGDMineta1_(1:MC,:);
        % -- 
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
        % 5
        SDAltGDFed = sum(SDAltGDFed,1)/MC;
        timeAltGDFed = sum(timeAltGDFed,1)/MC;
        % 7
        SDAltGDMineta1 = sum(SDAltGDMineta1,1)/MC;
        timeAltGDMineta1 = sum(timeAltGDMineta1,1)/MC;
        % 8
        SDAltGDMinCntrl = sum(SDAltGDMinCntrl,1)/MC;
        timeAltGDMinCntrl = sum(timeAltGDMinCntrl,1)/MC;
        %----------------------------------------------------------------------
        %------------ Downsample 
        % 2
        iterMax = numIter;
        SDAltMinCntrl = downsample(SDAltMinCntrl,1);
        timeAltMinCntrl = downsample(timeAltMinCntrl,1);
        iterAltMinCntrl = 1 : min(length(timeAltMinCntrl), iterMax+10);
        % 4
        SDAltGD = downsample(SDAltGD,4);
        timeAltGD = downsample(timeAltGD,4);
        iterAltGD = 4*(1 : min(length(timeAltGD), iterMax));
        % 7
        % 8
        SDAltGDMineta1 = downsample(SDAltGDMineta1,4);
        timeAltGDMineta1 = downsample(timeAltGDMineta1,4);   
        iterAltGDMineta1 = 4*(1 : min(length(timeAltGDMineta1), iterMax));
        %timeSD = 35;

        figure
        semilogy(iterAltGDMineta1,...
                SDAltGDMineta1(iterAltGDMineta1/4),'DisplayName', ...
               'AltGDMin', ...
               'LineWidth',1.35,'Marker','square','MarkerSize',9, 'Color',"#0000FF")
        hold on
        semilogy(iterAltMinCntrl,...
                SDAltMinCntrl(iterAltMinCntrl/1),...
               'DisplayName', ...
               'AltMin', ...
               'LineWidth',1.35,'Marker','x','MarkerSize',9,'Color',"#D95319")
        semilogy(iterAltGD,...
                 SDAltGD(iterAltGD/4),'DisplayName', ...
                'FactGD', ...
                'LineWidth',1.35,'Marker','diamond','MarkerSize',9,'Color',"#7E2F8E")
        grid on
        legend('Interpreter','Latex','Fontsize',9.25, 'Location','Northeast')
        yHndl = ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$','Interpreter','Latex','Fontsize',14.0);
        yHndl.Position(1) = yHndl.Position(1) + 1.25;
        yHndl.Position(2) = yHndl.Position(2) - 0.05*yHndl.Position(2);
        xlabel('Iterations','Fontsize',13,'Interpreter','Latex')
        title("n = " + n + ", q = " + q +...
              ", r = " + r + ", p = " + p + ", Workers =" + 10,...
               'Interpreter', 'Latex', 'FontSize',14)
        cores = feature('numCores');
        stringTitle = ['ErrAgnstIter_Cntrl_Wrkrs_',num2str(numWrkrs),'Max',num2str(cores),'_MC_',num2str(MC),...
                          '_n_',num2str(n),'_q_',num2str(q),'_r_',num2str(r),'_p_',num2str(p),...
                           'T_',num2str(T),'_id',num2str(randi(1e3,1))];
        exportgraphics(gcf,[stringTitle,'.pdf'])

end
