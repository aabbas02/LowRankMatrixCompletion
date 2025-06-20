function plotRsltsPhsTrnstn2d(numBlocks_, rVals, ...
                   prbAltGDMin, ...
                   prbAltMin,...
                   n,q,r,m,numBlocks,MC,same,T_LS,eta_c)
    curDir = pwd;
    cd ..
    cd figsMtrxSensing
    ID = num2str(randi(1e5));
    % ---- AltGDMin
    if any(any(prbAltGDMin))
        figure    
        %grid on        
        % Display using imagesc
        imagesc(prbAltGDMin);
        % x ticks
        xticks(1:length(numBlocks_));
        xticklabels(numBlocks_)
        ax = gca; % Get the current axis
        ax.XAxis.FontSize = 14; % Set the font size of the x-tick labels to 14
        % y ticks
        yticks(rVals)
        yticklabels(rVals)
        %set(gca,'yTickLabel','fontsize',12)    
        %gca.YAxis.FontSize = 25;
        ax = gca; % Get the current axis
        ax.YAxis.FontSize = 14; % Set the font size of the x-tick labels to 14
        % Set colormap to grayscale
        colormap(gray);
        % Ensure color scaling is [0, 1]
        caxis([0 1]);
        % Add colorbar for reference
        colorbar("FontSize",12);
        % Optional: remove axis ticks if you want a clean image display
        %axis off;
        %legend('Interpreter', 'Latex', 'Fontsize', 9);
        %ylabel("$\Pr[SD(U^(T),U^*)] < 10^{-10}$","FontSize",14,'Interpreter','Latex')
        xlabel('Number of blocks', 'FontSize',17, 'Interpreter','Latex')    
        ylabel('Rank of $X^*$', 'Interpreter','Latex','FontSize',17)
        title('AltGDMin. $\Pr [SD(U^{(T)},U^*) \leq 10^{-10}]$','Interpreter','Latex','FontSize',15)
        stringTitle = ['AltGDMin_PhaseTrnstn2d_MC_', num2str(MC), ...
                       '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), ... 
                       '_m_', num2str(m), '_numBlocks_',num2str(numBlocks), ...
                       '_same_',num2str(same),'_T_LS_', num2str(T_LS), '_eta_c_',num2str(eta_c), '_randID_',ID];
        
        savefig([stringTitle, '.fig']);
        exportgraphics(gca,[stringTitle,'.pdf'],"Resolution",300) 
    end
    % ---- AltMin
    if any(any(prbAltMin))
        figure    
        %grid on        
        % Display using imagesc
        imagesc(prbAltMin);
        % x ticks
        xticks(1:length(numBlocks_));
        xticklabels(numBlocks_)
        ax = gca; % Get the current axis
        ax.XAxis.FontSize = 14; % Set the font size of the x-tick labels to 14
        % y ticks
        yticks(rVals)
        yticklabels(rVals)
        %set(gca,'yTickLabel','fontsize',12)    
        %gca.YAxis.FontSize = 25;
        ax = gca; % Get the current axis
        ax.YAxis.FontSize = 14; % Set the font size of the x-tick labels to 14
        % Set colormap to grayscale
        colormap(gray);
        % Ensure color scaling is [0, 1]
        caxis([0 1]);
        % Add colorbar for reference
        colorbar("FontSize",12);
        % Optional: remove axis ticks if you want a clean image display
        %axis off;
        %legend('Interpreter', 'Latex', 'Fontsize', 9);
        %ylabel("$\Pr[SD(U^(T),U^*)] < 10^{-10}$","FontSize",14,'Interpreter','Latex')
        xlabel('Number of blocks', 'FontSize',17, 'Interpreter','Latex')    
        ylabel('Rank of $X^*$', 'Interpreter','Latex','FontSize',17)
        title('AltMin. $\Pr [SD(U^{(T)},U^*) \leq 10^{-10}]$','Interpreter','Latex','FontSize',15)
        stringTitle = ['AltMin_PhaseTrnstn2d_MC_', num2str(MC), ...
                       '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), ... 
                       '_m_', num2str(m), '_numBlocks_',num2str(numBlocks), ...
                       '_same_',num2str(same),'_T_LS_', num2str(T_LS),'_randID_',ID];
        
        savefig([stringTitle, '.fig']);
        exportgraphics(gca,[stringTitle,'.pdf'],"Resolution",300) 
    end
    cd (curDir)
end