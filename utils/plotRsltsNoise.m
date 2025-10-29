function plotRsltsNoise(time_AltGDMin, SDVals_AltGDMin, ...
                   time_UnPerm, SDVals_UnPerm, ...
                   time_AltMinExct,SDVals_AltMinExct,...
                   time_AltMin, SDVals_AltMin, ...
                   time_AltGDMinCllps, SDVals_AltGDMinCllps,...
                   time_AltMinExctCllps, SDVals_AltMinExctCllps,...
                   time_AltMinCllps, SDVals_AltMinCllps,...
                   n,q,r,m,numBlocks,MC,same,T_LS,eta_c,eta_L,noiseVars)
    curDir = pwd;
    cd ..
    cd figsMtrxSensing
    % AVERAGING
    %---
    SDVals_AltGDMin = squeeze(sum(SDVals_AltGDMin,2)/MC);
    SDVals_UnPerm = squeeze(sum(SDVals_UnPerm,2)/MC);
    SDVals_AltMinExct = squeeze(sum(SDVals_AltMinExct,2)/MC);   
    SDVals_AltMin = squeeze(sum(SDVals_AltMin,2)/MC);
    %---
    SDVals_AltGDMinCllps = squeeze(sum(SDVals_AltGDMinCllps,2)/MC);
    SDVals_AltMinExctCllps = squeeze(sum(SDVals_AltMinExctCllps,2)/MC);   
    SDVals_AltMinCllps = squeeze(sum(SDVals_AltMinCllps,2)/MC);
    %------------------------------------
    time_UnPerm  = squeeze(sum(time_UnPerm,2)/MC);
    time_AltGDMin = squeeze(sum(time_AltGDMin,2)/MC);
    time_AltMinExct = squeeze(sum(time_AltMinExct,2)/MC);
    time_AltMin = squeeze(sum(time_AltMin,2)/MC);
    %---
    time_AltGDMinCllps = squeeze(sum(time_AltGDMinCllps,2)/MC);
    time_AltMinExctCllps = squeeze(sum(time_AltMinExctCllps,2)/MC);
    time_AltMinCllps = squeeze(sum(time_AltMinCllps,2)/MC);    
    %--------------------------------------
    % SD vs iter figure
    %--------------------------------------
    %{
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
        for i = 1 : len(noiseVars)
            noiseVar = noiseVars(i);
            semilogy(SDVals_AltMinExct, ...
                'DisplayName', 'AltMin (Exact)' + num2str(noiseVar), 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 7);
        end
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
                   '_noiseVar_',num2str(noiseVar),...
                   'Block_Size_',num2str(m/numBlocks), '_same_',num2str(same), '_T_LS_', num2str(T_LS)];
    savefig([stringTitle, '.fig']);   
    %}
    %--------------------
    % SD VS TIME FIGURE
    %--------------------
    figure
    if ~(all (SDVals_AltGDMin == 0))
        for i = 1 : length(noiseVars)
            semilogy(time_AltGDMin(i,:), SDVals_AltGDMin(i,:), ...
                'DisplayName', ['AltGDMin $\sigma^2 = $' num2str(noiseVars(i))], 'LineWidth', 1.45, 'Marker', 'diamond', 'MarkerSize', 7);
            hold on;
        end
    end
    if ~(all (SDVals_UnPerm == 0))
        for i = 1 : length(noiseVars)
            semilogy(time_UnPerm(i,:), SDVals_UnPerm(i,:), ...
                'DisplayName', ['Unpermuted $\sigma^2 = $'  num2str(noiseVars(i))], 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
            hold on
        end
    end
    if ~(all (SDVals_AltMinExct == 0))
        for i = 1 : length(noiseVars)
            semilogy(time_AltMinExct(i,:), SDVals_AltMinExct(i,:), ...
                'DisplayName', ['AltMin (Exact) $\sigma^2 = $'  num2str(noiseVars(i))], 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);    
            hold on
        end
    end
    if ~(all (SDVals_AltMin == 0))
        for i = 1: length(noiseVars)
            semilogy(time_AltMin(i,:), SDVals_AltMin(i,:),...
            'DisplayName', ['AltMin (GD) $\sigma^2 = $'  num2str(noiseVars(i))], 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7)
            hold on
        end
    end
    % Collapsed
    if ~(all (SDVals_AltGDMinCllps == 0))
        for i = 1 : length(noiseVars)
            semilogy(time_AltGDMinCllps(i,:), SDVals_AltGDMinCllps(i,:), ...
                'DisplayName', [' AltGDMin - Cllps $\sigma^2 = $' num2str(noiseVars(i))], 'LineWidth', 1.45, 'Marker', 'diamond', 'MarkerSize', 7);
            hold on;
        end
    end
    if ~(all (SDVals_UnPerm == 0))
        for i = 1 : length(noiseVars)
            semilogy(time_UnPerm(i,:), SDVals_UnPerm(i,:), ...
                'DisplayName', ['Unpermuted $\sigma^2 = $' num2str(noiseVars(i))], 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
            hold on
        end
    end
    if ~(all (SDVals_AltMinExctCllps == 0))
        for i = 1 : length(noiseVars)
            semilogy(time_AltMinExctCllps(i,:), SDVals_AltMinExctCllps(i,:), ...
                'DisplayName', ['AltMin (Exact) - Cllps $\sigma^2 = $'  num2str(noiseVars(i))], 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);    
            hold on
        end
    end
    if ~(all (SDVals_AltMinCllps == 0))
        for i = 1 : length(noiseVars)
            semilogy(time_AltMinCllps(i,:), SDVals_AltMinCllps(i,:),...
                 'DisplayName', ['AltMin (GD) - Cllps $\sigma^2 = $' num2str(noiseVars(i))], 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7)
                  hold on
        end
    end
    legend('show', 'Interpreter', 'latex');
    grid on    
    ax = gca();
    ax.FontSize = 12;
    title("n = " + n + ", q = " + q +...
          ", r = " + r + ", m = " + m + ", Block Size = " ...
          + m/numBlocks +  ", MC = " + MC + ", $\eta_c$ = " + eta_c, ...
         'Interpreter', 'Latex', 'FontSize',11)
    
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