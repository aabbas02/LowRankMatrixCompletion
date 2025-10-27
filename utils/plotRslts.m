function plotRslts(time_AltGDMin, SDVals_AltGDMin, ...
                   time_UnPerm, SDVals_UnPerm, ...
                   time_AltMinExct,SDVals_AltMinExct,...
                   time_AltMin, SDVals_AltMin, ...
                   time_AltGDMinCllps, SDVals_AltGDMinCllps,...
                   time_AltMinExctCllps, SDVals_AltMinExctCllps,...
                   time_AltMinCllps, SDVals_AltMinCllps,...
                   n,q,r,m,numBlocks,MC,same,T_LS,eta_c,eta_L)
    curDir = pwd;
    cd ..
    cd figsMtrxSensing
    % AVERAGING
    %---
    SDVals_AltGDMin = sum(SDVals_AltGDMin,1)/MC;
    SDVals_UnPerm = sum(SDVals_UnPerm,1)/MC;
    SDVals_AltMinExct = sum(SDVals_AltMinExct,1)/MC;   
    SDVals_AltMin = sum(SDVals_AltMin,1)/MC;
    %---
    SDVals_AltGDMinCllps = sum(SDVals_AltGDMinCllps,1)/MC;
    SDVals_AltMinExctCllps = sum(SDVals_AltMinExctCllps,1)/MC;   
    SDVals_AltMinCllps = sum(SDVals_AltMinCllps,1)/MC;
    %------------------------------------
    time_UnPerm  = sum(time_UnPerm,1)/MC;
    time_AltGDMin = sum(time_AltGDMin,1)/MC;
    time_AltMinExct = sum(time_AltMinExct,1)/MC;
    time_AltMin = sum(time_AltMin,1)/MC;
    %---
    time_AltGDMinCllps = sum(time_AltGDMinCllps,1)/MC;
    time_AltMinExctCllps = sum(time_AltMinExctCllps,1)/MC;
    time_AltMinCllps = sum(time_AltMinCllps,1)/MC;    
    %--------------------------------------
    % SD vs iter figure
    %--------------------------------------
    figure;
    if ~(all (SDVals_AltGDMin == 0))
        semilogy(SDVals_AltGDMin, ...
            'DisplayName', 'AltGDMin', 'LineWidth', 1.45, 'Marker', 'diamond', 'MarkerSize', 7);
        hold on;
    end
    if ~(all (SDVals_UnPerm == 0))
        semilogy(SDVals_UnPerm, ...
            'DisplayName', 'Unpermuted', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
        hold on;        
    end
    if ~(all (SDVals_AltMinExct == 0))
        semilogy(SDVals_AltMinExct, ...
            'DisplayName', 'AltMin (Exact)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 7);
        hold on
    end
    if ~(all (SDVals_AltMin == 0))    
        semilogy(SDVals_AltMin,...
             'DisplayName', 'AltMin (GD)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7)
        hold on
    end
    %--------------- Cllps
    if ~(all (SDVals_AltGDMinCllps == 0))
        semilogy(SDVals_AltGDMinCllps, ...
            'DisplayName', 'AltGDMin - Cllps', 'LineWidth', 1.45, 'Marker', 'diamond', 'MarkerSize', 7);
        hold on;
    end
    if ~(all (SDVals_UnPerm == 0))
        semilogy(SDVals_UnPerm, ...
            'DisplayName', 'Unpermuted', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
        hold on;        
    end
    if ~(all (SDVals_AltMinExctCllps == 0))
        semilogy(SDVals_AltMinExctCllps, ...
            'DisplayName', 'AltMin (Exact) - Cllps', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 7);
        hold on
    end
    if ~(all (SDVals_AltMinCllps == 0))    
        semilogy(SDVals_AltMinCllps,...
             'DisplayName', 'AltMin (GD) - Cllps', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7)
        hold on
    end 
    ax = gca();
    ax.FontSize = 12;
    grid on    
    title("n = " + n + ", q = " + q +...
          ", r = " + r + ", m = " + m + ", Block Size = " + m/numBlocks +  ", MC = " + MC, ...
           'Interpreter', 'Latex', 'FontSize',12)
    
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel("$SD(U^{(t)},U^*)$","FontSize",14,'Interpreter','Latex')
    xlabel('Iterations (t)', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['Iter_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), ...
                   '_r_', num2str(r), '_m_', num2str(m), ...
                   'Block_Size_',num2str(m/numBlocks), '_same_',num2str(same), '_T_LS_', num2str(T_LS)];
    savefig([stringTitle, '.fig']);   
    %--------------------
    % SD VS TIME FIGURE
    %--------------------
    figure
    if ~(all (SDVals_AltGDMin == 0))
    semilogy(time_AltGDMin, SDVals_AltGDMin, ...
        'DisplayName', ' AltGDMin', 'LineWidth', 1.45, 'Marker', 'diamond', 'MarkerSize', 7);
    hold on;
    end
    if ~(all (SDVals_UnPerm == 0))
        semilogy(time_UnPerm, SDVals_UnPerm, ...
            'DisplayName', 'Unpermuted', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
        hold on
    end
    if ~(all (SDVals_AltMinExct == 0))
        semilogy(time_AltMinExct, SDVals_AltMinExct, ...
            'DisplayName', 'AltMin (Exact)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);    
        hold on
    end
    if ~(all (SDVals_AltMin == 0))
    semilogy(time_AltMin, SDVals_AltMin,...
         'DisplayName', 'AltMin (GD)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7)
          hold on
    end
    % Collapsed
    if ~(all (SDVals_AltGDMinCllps == 0))
    semilogy(time_AltGDMinCllps, SDVals_AltGDMinCllps, ...
        'DisplayName', ' AltGDMin - Cllps', 'LineWidth', 1.45, 'Marker', 'diamond', 'MarkerSize', 7);
    hold on;
    end
    if ~(all (SDVals_UnPerm == 0))
        semilogy(time_UnPerm, SDVals_UnPerm, ...
            'DisplayName', 'Unpermuted', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
        hold on
    end
    if ~(all (SDVals_AltMinExctCllps == 0))
        semilogy(time_AltMinExctCllps, SDVals_AltMinExctCllps, ...
            'DisplayName', 'AltMin (Exact) - Cllps', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);    
        hold on
    end
    if ~(all (SDVals_AltMinCllps == 0))
    semilogy(time_AltMinCllps, SDVals_AltMinCllps,...
         'DisplayName', 'AltMin (GD) - Cllps', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7)
          hold on
    end
    grid on    
    ax = gca();
    ax.FontSize = 12;
    title("n = " + n + ", q = " + q +...
          ", r = " + r + ", m = " + m + ", Block Size = " ...
          + m/numBlocks +  ", MC = " + MC + ", $\eta_c$ = " + eta_c, ...
           'Interpreter', 'Latex', 'FontSize',12)
    
    legend('Interpreter', 'Latex', 'Fontsize', 11);
    ylabel("$SD(U^{(t)},U^*)$","FontSize",14,'Interpreter','Latex')
    xlabel('t (Seconds)', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['Time_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), ... 
                   '_m_', num2str(m), '_blockSize_',num2str(m/numBlocks), ...
                   '_same_',num2str(same),'_T_LS_', num2str(T_LS), '_eta_c_',num2str(eta_c),...
                   '_eta_L_',num2str(eta_L)];
    savefig([stringTitle, '.fig']);    
    exportgraphics(gca,[stringTitle,'.pdf'],"Resolution",300)
    cd (curDir)
end