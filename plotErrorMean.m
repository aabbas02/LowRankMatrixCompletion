function plotErrorMean(error_mean_MC, n, q, r, p, L, L_byz, MC, C1,attck)
    % Create a figure
    figure;
    % Plot the errors
    %semilogy(error_mean_MC(1,:), 'Color', "#7E2F8E",'LineWidth',1.25,'DisplayName',"Sum (No Attack)");
    len = length(error_mean_MC(2,:));
    if len > 40
        fctr = len/40;
        T_ = 1:len;
        T_ = downsample(T_,floor(fctr));
        semilogy(T_, error_mean_MC(2,T_), 'Color', '#0072BD', 'LineWidth', 2.5, 'DisplayName', "Elmntwise Mdn.", 'Marker', 'x', 'MarkerSize', 11);
        hold on
        semilogy(T_, error_mean_MC(4,T_), 'Color', '#EDB120', 'LineWidth', 2.5, 'DisplayName', "Krum",'Marker', 'square', 'MarkerSize', 11);
    else
        semilogy(error_mean_MC(2,:), 'Color', '#0072BD', 'LineWidth', 2.5, 'DisplayName', "Elmntwise Median", 'Marker', 'x', 'MarkerSize', 11);
        hold on;
        %semilogy(error_mean_MC(3,:), 'Color', '#D95319','LineWidth',1.25, 'DisplayName',"Goemetric Mdn.");
        semilogy(error_mean_MC(4,:), 'Color', '#EDB120', 'LineWidth', 2.5, 'DisplayName', "Krum",'Marker', 'square', 'MarkerSize', 11);
    end
    %hold on;
    %semilogy(error_mean_MC(3,:), 'Color', '#D95319','LineWidth',1.25, 'DisplayName',"Goemetric Mdn.");
    legend('Interpreter', 'Latex', 'Fontsize', 11,'Location','southwest');
    % Generate a random ID
    ID = randi(10000);
    grid("on")
    % Set the title with dynamic parameters
    xlabel('Iterations $t$', 'Interpreter', 'Latex','FontSize',14);
    ylabel('$\mathrm{SD}(U^{(t)},U^*)$', 'Interpreter', 'Latex','FontSize',14);
    title(["$n = $ " + n + ", $ q = $" + q + ...
           "$, r = $ " + r + ", $ p = $ " + p + ...
           "$, L = $ " + L + ", $L_{byz} = $ " + L_byz + ...
           ", $MC = $" + num2str(MC) +  ",  Attck = " + num2str(attck) + "."], ...
          'Interpreter', 'Latex', 'FontSize', 12)
    % Create a string for the filename
    stringTitle = ['New_MC_', num2str(MC), ...
                   '_n_', num2str(n), '_q_', num2str(q), '_r_', ...
                   num2str(r), '_L_', num2str(L), 'p', num2str(p), ...
                   '_Lbyz_', num2str(L_byz), 'C1_', num2str(C1), 'attck_', num2str(attck) ...
                   '_ID', num2str(ID)];
    % Save the figure
    savefig([stringTitle, '.fig']);
    save("vars"+ stringTitle + ".mat","error_mean_MC","n","q","r","p","L","L_byz","MC","C1","attck")
end