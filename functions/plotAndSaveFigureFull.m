function plotAndSaveFigureFull(timeAltMinCntrl, SDAltMinCntrl, timeAltMinParfor, SDAltMinParfor, ...
                               timeAltMinPrvt, SDAltMinPrvt, timeAltGD, SDAltGD, timeAltGDFed, SDAltGDFed, ...
                               timeAltGDMineta1, SDAltGDMineta1, timeAltGDMinCntrl, SDAltGDMinCntrl,  ...
                                n, q, r, p, numWrkrs, MC, T,Tsvd)

    idx1 = find(SDAltMinParfor < 10^-14,1);
    tidx1 = timeAltMinParfor(idx1);

    idx2 = find(SDAltMinCntrl < 10^-14,1);
    tidx2 = timeAltMinCntrl(idx2);

    idx3 = find(SDAltMinPrvt < 10^-14,1);
    tidx3 = timeAltMinPrvt(idx3);

    idx4 =  find(SDAltGD < 10^-14,1);
    tidx4 = timeAltGD(idx4);

    idx5 = find(SDAltGDFed  < 10^-15,1);
    tidx5 = timeAltGDFed(idx5);

    %idx6 = find(SDAltGDMin < 10^-14,1);
    %tidx6 =  timeAltGDMin(idx6);

    idx7 = find(SDAltGDMineta1 < 10^-14,1);
    tidx7 =  timeAltGDMineta1(idx7);

    idx8 = find(SDAltGDMinCntrl < 10^-14,1);
    tidx8 =  timeAltGDMinCntrl(idx8);

    %timeSD = max([tidx1,tidx2,tidx3,tidx4,tidx5,tidx6,tidx7,tidx8]);
    timeSD = max([tidx1,tidx2,tidx3,tidx4,tidx5,tidx7,tidx8]);
    %timeSD = max([tidx1,tidx2,tidx4,tidx6,tidx7,tidx8]);
    %timeSD = max([tidx1,tidx7]);


    figure;
    semilogy(timeAltMinCntrl(timeAltMinCntrl <= timeSD), ...
        SDAltMinCntrl(timeAltMinCntrl <= timeSD), ...
        'DisplayName', 'AltMin(Cntrl.)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 5);
    hold on;
    semilogy(timeAltMinParfor(timeAltMinParfor <= timeSD), ...
        SDAltMinParfor(timeAltMinParfor <= timeSD), ...
        'DisplayName', 'AltMin(Fed./NotPrvt.)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 5);
    semilogy(timeAltMinPrvt(timeAltMinPrvt <= timeSD), ...
        SDAltMinPrvt(timeAltMinPrvt <= timeSD), ...
        'DisplayName', 'AltMin(Fed./Prvt.)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 5);
    semilogy(timeAltGD(timeAltGD <= timeSD), ...
        SDAltGD(timeAltGD <= timeSD), ...
        'DisplayName', 'FactGD (Cntrl.)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    semilogy(timeAltGDFed(timeAltGDFed <= timeSD), ...
        SDAltGDFed(timeAltGDFed <= timeSD), ...
        'DisplayName', 'FactGD (Fed.)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 5);
    semilogy(timeAltGDMineta1(timeAltGDMineta1 <= timeSD), ...
        SDAltGDMineta1(timeAltGDMineta1 <= timeSD), ...
        'DisplayName', 'AltGDMin(Fed.)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 5);
    semilogy(timeAltGDMinCntrl(timeAltGDMinCntrl <= timeSD), ...
        SDAltGDMinCntrl(timeAltGDMinCntrl <= timeSD), ...
        'DisplayName', 'AltGDMin(Cntrl.)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 5);
    
    grid on;
    
    legend('Interpreter', 'Latex', 'Fontsize', 9);
    ylabel('$\mathrm{SD}(\mathbf{U}^{(t)}, \mathbf{U}^*)$', 'Interpreter', 'Latex', 'Fontsize', 15);
    xlabel('Time/seconds', 'FontSize', 11);
    
    title("n = " + n + ", q = " + q +...
        ", r = " + r + ", p = " + p + '.', ...
        'Interpreter', 'Latex', 'FontSize', 14);
    
    cores = feature('numCores');
    %stringTitle = ['Wrkrs', num2str(numWrkrs),'_MC_', num2str(MC), ...
    %    '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_p_', num2str(p), ...
    %    'T_', num2str(T),'_Tsvd',num2str(Tsvd),'_id', num2str(randi(1e3, 1))];
    
    %savefig([stringTitle, '.fig']);
end
