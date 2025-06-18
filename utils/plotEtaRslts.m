function plotEtaRslts(eta1,time_AltGDMin1, SDVals_AltGDMin1, ...
                   eta2, time_AltGDMin2, SDVals_AltGDMin2,...
                   eta3, time_AltGDMin3, SDVals_AltGDMin3,...
                   time_AltMin, SDVals_AltMin,...
                   n,q,r,m,numBlocks,MC,same,T_LS,eta_c)
    curDir = pwd;
    cd ..
    cd figsMtrxSensing
    %---
    SDVals_AltGDMin1 = sum(SDVals_AltGDMin1,1)/MC;
    SDVals_AltGDMin2 = sum(SDVals_AltGDMin2,1)/MC;
    SDVals_AltGDMin3 = sum(SDVals_AltGDMin3,1)/MC;
    SDVals_AltMin = sum(SDVals_AltMin,1)/MC;
    %---
    time_AltGDMin1 = sum(time_AltGDMin1,1)/MC;
    time_AltGDMin2 = sum(time_AltGDMin2,1)/MC;
    time_AltGDMin3 = sum(time_AltGDMin3,1)/MC;
    time_AltMin = sum(time_AltMin,1)/MC;
    figure
    if ~(all (SDVals_AltGDMin1 == 0))
    semilogy(time_AltGDMin1, SDVals_AltGDMin1, ...
        'DisplayName', ['$\eta_c = $ ', num2str(eta1)], 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 5);
    hold on;
    end
    if ~(all (SDVals_AltGDMin2 == 0))
        semilogy(time_AltGDMin2, SDVals_AltGDMin2, ...
            'DisplayName', ['$\eta_c = $ ', num2str(eta2)], 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 5);
        hold on
    end
    if ~(all (SDVals_AltGDMin3 == 0))
        semilogy(time_AltGDMin3, SDVals_AltGDMin3, ...
            'DisplayName', ['$\eta_c = $ ', num2str(eta3)], 'LineWidth', 1.45, 'Marker', '+', 'MarkerSize', 5);    
        hold on
    end
    if ~(all (SDVals_AltMin == 0))
        semilogy(time_AltMin, SDVals_AltMin, ...
            'DisplayName', 'AltMin (GD)', 'LineWidth', 1.45, 'Marker', 'diamond', 'MarkerSize', 5);    
        hold on
    end    
    grid on    
    title("n = " + n + ", q = " + q +...
          ", r = " + r + ", m = " + m + ", num Blocks = " ...
          + numBlocks +  ", MC = " + MC + ", $\eta_c$ = " + eta_c, ...
           'Interpreter', 'Latex', 'FontSize',11)
    
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel("$SD(U^{(t)},U^*)$","FontSize",14,'Interpreter','Latex')
    xlabel('t (Seconds)', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['Time_Eta_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), ... 
                   '_m_', num2str(m), '_numBlocks_',num2str(numBlocks), ...
                   '_same_',num2str(same),'_T_LS_', num2str(T_LS), '_eta_c_',num2str(eta_c)];
    savefig([stringTitle, '.fig']);    
    exportgraphics(gca,[stringTitle,'.pdf'],"Resolution",300)
    cd (curDir)
end