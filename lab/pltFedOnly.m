function pltFedOnly(timeAltMinParfor, SDAltMinParfor, ...
                    timeAltMinPrvt, SDAltMinPrvt, timeAltGDFedHalf, SDAltGDFedHalf, ...
                    timeAltGDFed,SDAltGDFed, timeAltGDFed1,SDAltGDFed1,...
                    timeAltGDMineta1, SDAltGDMineta1,...
                    n, q, r, p, numWrkrs, MC, T,Tsvd)

    idx1 = find(SDAltMinParfor < 10^-14,1);
    tidx1 = timeAltMinParfor(idx1);
    %
    idx3 = find(SDAltMinPrvt < 10^-14,1);
    tidx3 = timeAltMinPrvt(idx3);
    %
    idx4 = find(SDAltGDFedHalf  < 10^-15,1);
    tidx4 = timeAltGDFedHalf(idx4);
    %
    idx5 = find(SDAltGDFed  < 10^-15,1);
    tidx5 = timeAltGDFed(idx5);
    %
    idx6 = find(SDAltGDFed1  < 10^-15,1);
    tidx6 = timeAltGDFed1(idx6);
    %
    idx7 = find(SDAltGDMineta1 < 10^-14,1);
    tidx7 =  timeAltGDMineta1(idx7);
    %
    timeSD = max([tidx1,tidx3,tidx4,tidx5,tidx6,tidx7]);
    %
    figure;
    semilogy(timeAltMinParfor(timeAltMinParfor <= timeSD), ...
        SDAltMinParfor(timeAltMinParfor <= timeSD), ...
        'DisplayName', 'AltMin(Fed./NotPrvt.)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 5);
    hold on;
    semilogy(timeAltMinPrvt(timeAltMinPrvt <= timeSD), ...
        SDAltMinPrvt(timeAltMinPrvt <= timeSD), ...
        'DisplayName', 'AltMin(Fed./Prvt.)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 5);
    %
    semilogy(timeAltGDFedHalf(timeAltGDFedHalf <= timeSD), ...
        SDAltGDFedHalf(timeAltGDFedHalf <= timeSD), ...
        'DisplayName', 'FactGD (Fed. $c = 0.5$)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    %
    semilogy(timeAltGDFed(timeAltGDFed <= timeSD), ...
        SDAltGDFed(timeAltGDFed <= timeSD), ...
        'DisplayName', 'FactGD (Fed. $c = 0.75$)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    %
    semilogy(timeAltGDFed1(timeAltGDFed1 <= timeSD), ...
        SDAltGDFed(timeAltGDFed1 <= timeSD), ...
        'DisplayName', 'FactGD (Fed. $c = 1.00$)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    % 
    semilogy(timeAltGDMineta1(timeAltGDMineta1 <= timeSD), ...
        SDAltGDMineta1(timeAltGDMineta1 <= timeSD), ...
        'DisplayName', 'AltGDMin(Fed.)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 5);
    % 
    grid on;
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$', 'Interpreter', 'Latex', 'Fontsize', 15);
    xlabel('Time/seconds', 'FontSize', 11);
    
    title("n = " + n + ", q = " + q +...
        ", r = " + r + ", p = " + p + '.', ...
        'Interpreter', 'Latex', 'FontSize', 14);
    %
    %cores = feature('numCores');
    %stringTitle = ['FedOnly_', num2str(numWrkrs),'_MC_', num2str(MC), ...
    %    '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_p_', num2str(p), ...
    %    'T_', num2str(T),'_Tsvd',num2str(Tsvd),'_id', num2str(randi(1e3, 1))];
    
    %savefig([stringTitle, '.fig']);
end
