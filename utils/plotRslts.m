function plotRslts(time_AltGDMin, SDVals_AltGDMin, ...
                   time_UnPerm, SDVals_UnPerm, ...
                   time_AltMin, SDVals_AltMin, ...
                   n,q,r,m,numBlocks,MC,same)

    SDVals_UnPerm = sum(SDVals_UnPerm,1)/MC;
    SDVals_AltGDMin = sum(SDVals_AltGDMin,1)/MC;
    SDVals_AltMin = sum(SDVals_AltMin,1)/MC;
    %---
    time_UnPerm  = sum(time_UnPerm,1)/MC;
    time_AltGDMin = sum(time_AltGDMin,1)/MC;
    time_AltMin = sum(time_AltMin,1)/MC;
    % SD vs iter figure
    figure;
    
    semilogy(SDVals_AltGDMin, ...
        'DisplayName', 'SDVals (s-Local)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    hold on;
    %semilogy(SDVals_Perm, ...
    %    'DisplayName', 'SDVals (Naive)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 7);
    %
    semilogy(SDVals_UnPerm, ...
        'DisplayName', 'SDVals (Unpermuted)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
    semilogy(SDVals_AltMin,...
         'DisplayName', 'SDVals (AltMin)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7)
    grid on
    
    
    title("n = " + n + ", q = " + q +...
          ", r = " + r + ", m = " + m + ", num Blocks = " + numBlocks +  ", MC = " + MC + ", same = " + same, ...
           'Interpreter', 'Latex', 'FontSize',14)
    
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel("$SD(U,U^*)$","FontSize",14,'Interpreter','Latex')
    xlabel('Iterations (t)', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['Iter_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_m_', num2str(m), '_numBlocks_',num2str(numBlocks), '_same_',num2str(same)];
    
    savefig([stringTitle, '.fig']);
    figure
    SDVals_UnPerm = sum(SDVals_UnPerm,1)/MC;
    SDVals_AltMin = sum(SDVals_AltMin,1)/MC;
    
    semilogy(time_AltGDMin, SDVals_AltGDMin, ...
        'DisplayName', 'SDVals (s-Local)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    hold on;
    %semilogy(SDVals_Perm, ...
    %    'DisplayName', 'SDVals (Naive)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 7);
    %
    semilogy(time_UnPerm, SDVals_UnPerm, ...
        'DisplayName', 'SDVals (Unpermuted)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
    semilogy(time_AltMin, SDVals_AltMin,...
         'DisplayName', 'SDVals (AltMin)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7)
    grid on
    
    
    title("n = " + n + ", q = " + q +...
          ", r = " + r + ", m = " + m + ", num Blocks = " + numBlocks +  ", MC = " + MC + ", same = " + same, ...
           'Interpreter', 'Latex', 'FontSize',14)
    
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel("$SD(U,U^*)$","FontSize",14,'Interpreter','Latex')
    xlabel('Seconds (Time)', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['Time_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_m_', num2str(m), '_numBlocks_',num2str(numBlocks), '_same_',num2str(same)];
    
    savefig([stringTitle, '.fig']);

end