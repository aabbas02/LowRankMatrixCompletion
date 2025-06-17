function plotRsltsPhsTrnstn2d(numBlocks_, rVals, ...
                   prbAltGDMin, ...
                   prbAltMin,...
                   n,q,r,m,numBlocks,MC,same,T_LS,eta_c)
    curDir = pwd;
    cd ..
    cd figsMtrxSensing
    figure
    %grid on        

    % Display using imagesc
    imagesc(prbAltGDMin);

    % x ticks
    xticks(numBlocks_);

    % y ticks
    yticks(rVals)
    
    % Set colormap to grayscale
    colormap(gray);
    
    % Ensure color scaling is [0, 1]
    caxis([0 1]);
    
    % Add colorbar for reference
    colorbar;
    
    % Optional: remove axis ticks if you want a clean image display
    axis off;
    title("AltGDMin. n = " + n + ", q = " + q +...
          ", r = " + r + ", m = " + m + ", MC = " + MC + ", $\eta_c = $ " + eta_c, ...
           'Interpreter', 'Latex', 'FontSize',11)
    
    %legend('Interpreter', 'Latex', 'Fontsize', 9);
    %ylabel("$\Pr[SD(U^(T),U^*)] < 10^{-10}$","FontSize",14,'Interpreter','Latex')
    xlabel('Number of blocks', 'FontSize',14, 'Interpreter','Latex')    
    stringTitle = ['PhaseTrnstn2d_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), ... 
                   '_m_', num2str(m), '_numBlocks_',num2str(numBlocks), ...
                   '_same_',num2str(same),'_T_LS_', num2str(T_LS), '_eta_c_',num2str(eta_c)];
    
    savefig([stringTitle, '.fig']);
    exportgraphics(gca,[stringTitle,'.pdf'],"Resolution",300)
    cd (curDir)
end