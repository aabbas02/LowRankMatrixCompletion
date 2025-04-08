function plotRsltsLRMC(SDU0, SDU0Cllps, SDU0Cllps_, SDU0Perm,n,q,r,p,numBlocks_,MC,same,fill,real,T)
    figure;
    hold on
    if real == 0
    plot(numBlocks_,SDU0, ...
        'DisplayName', 'SDVals Unperm', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    end
    %plot(numBlocks_,SDU0Cllps, ...
    %    'DisplayName', 'SDVals Cllps', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    plot(numBlocks_,SDU0Cllps_, ...
        'DisplayName', 'SDVals Cllps Mean Only', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    %plot(numBlocks_,SDU0Perm, ...
    %    'DisplayName', 'SDVals Naive', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
    grid("on")
    xticks(numBlocks_);
    %-------------------------------
    title("LRMC. n = " + n + ", q = " + q +...
          ", r = " + r + ", p = " + p  +  ", MC = " + MC + ", same = " + same + ", fill = " + fill + ", T = " + T,...
           'Interpreter', 'Latex', 'FontSize',11)
    %--------------------------------
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    if real == 0
        ylabel("$SD(U^{(T)},U^*)$","FontSize",14,'Interpreter','Latex')
    else
        ylabel("Initialization Error", "FontSize",11,"Interpreter","Latex")
    end
    xlabel('number of blocks', 'FontSize',14, 'Interpreter','Latex')
    stringTitle = ['LRMC_Init_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_p_', num2str(p),'_numBlocks_', num2str(max(numBlocks_)), '_same_',num2str(same),'_fill_',num2str(fill), '_T_',num2str(T)];
    
    savefig([stringTitle, '.fig']);
end
%---
