function plotRsltsPhsTrnstn(numBlocks_, ...
                   prbAltGDMin, ...
                   prbAltMin,...
                   n,q,r,m,numBlocks,MC,same,T_LS,eta_c)
    curDir = pwd;
    cd ..
    cd figsMtrxSensing
    % SD vs iter figure
    figure;
    %{
    if ~(all (SDVals_AltGDMin == 0))
        semilogy(SDVals_AltGDMin, ...
            'DisplayName', 'AltGDMin', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 1);
        hold on;
    end
    if ~(all (SDVals_UnPerm == 0))
        semilogy(SDVals_UnPerm, ...
            'DisplayName', 'Unpermuted', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 1);
        hold on;        
    end
    if ~(all (SDVals_AltMinExct == 0))
        semilogy(SDVals_AltMinExct, ...
            'DisplayName', 'AltMin (Exact)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 1);
        hold on
    end
    if ~(all (SDVals_AltMin == 0))    
        semilogy(SDVals_AltMin,...
             'DisplayName', 'AltMin (GD)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 1)
        hold on
    end
    %}
    plot(numBlocks_, prbAltGDMin, 'DisplayName', 'AltGDMin', 'LineWidth', 1.05, 'Marker', 'o', 'MarkerSize', 3)
    hold on
    plot(numBlocks_, prbAltMin, 'DisplayName', 'AltMin', 'LineWidth', 1.05, 'Marker', 'square', 'MarkerSize', 3)
    grid on        
    title("n = " + n + ", q = " + q +...
          ", r = " + r + ", m = " + m + ", MC = " + MC + ", $\eta_c = $ " + eta_c, ...
           'Interpreter', 'Latex', 'FontSize',11)
    
    legend('Interpreter', 'Latex', ['' ...
        'Fontsize'], 9);
    ylabel("$\Pr[SD(U^(T),U^*)] < 10^{-10}$","FontSize",14,'Interpreter','Latex')
    xlabel('Number of blocks', 'FontSize',14, 'Interpreter','Latex')    
    stringTitle = ['PhaseTrnstn_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), ... 
                   '_m_', num2str(m), '_numBlocks_',num2str(numBlocks), ...
                   '_same_',num2str(same),'_T_LS_', num2str(T_LS), '_eta_c_',num2str(eta_c)];
    
    savefig([stringTitle, '.fig']);
    exportgraphics(gca,[stringTitle,'.pdf'],"Resolution",300)
    cd (curDir)
end