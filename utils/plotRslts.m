function plotRslts(SDVals_sLcl, SDVals_Perm, SDVals_UnPerm,n,q,r,m,numBlocks,MC,same)
    figure;
    SDVals_UnPerm = sum(SDVals_UnPerm,1)/MC;
    SDVals_sLcl = sum(SDVals_sLcl,1)/MC;
    semilogy(SDVals_sLcl, ...
        'DisplayName', 'SDVals (s-Local)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    hold on;
    %semilogy(SDVals_Perm, ...
    %    'DisplayName', 'SDVals (Naive)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 7);
    %
    semilogy(SDVals_UnPerm, ...
        'DisplayName', 'SDVals (Unpermuted)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
    grid on
    
    
    title("n = " + n + ", q = " + q +...
          ", r = " + r + ", m = " + m + ", num Blocks = " + numBlocks +  ", MC = " + MC + ", same = " + same, ...
           'Interpreter', 'Latex', 'FontSize',14)
    
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel("$SD(U,U^*)$","FontSize",14,'Interpreter','Latex')
    xlabel('Iterations (t)', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['MtrxSnsngPerm_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_m_', num2str(m), '_numBlocks_',num2str(numBlocks), '_same_',num2str(same)];
    
    savefig([stringTitle, '.fig']);
end