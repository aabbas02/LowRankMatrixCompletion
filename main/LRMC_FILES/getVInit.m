function V_init = getVInit(r, q, r_All, rowIdx, U, Xcol)
    % Initialize V_init matrix
    V_init = zeros(r, q);
    % Loop over each column of the V_init matrix
    for k = 1 : q
        % Initialize Ucllpsk and XcolCllpsk
        Ucllpsk = zeros(length(r_All{k}), r);
        XcolCllpsk = zeros(length(r_All{k}), 1);
        start = 1;        
        % Loop over each block in r_All{k}
        for i = 1 : length(r_All{k})
            blkSize = r_All{k}(i);
            rowsBlock = rowIdx{k}(start : start + blkSize - 1);            
            % Update Ucllpsk and XcolCllpsk
            Ucllpsk(i, :) = sum(U(rowsBlock, :));
            XcolCllpsk(i) = sum(Xcol{k}(start : start + blkSize - 1));
            start = start + blkSize;
        end
        % Least-squares update of V_init
        V_init(:, k) = pinv(Ucllpsk) * XcolCllpsk;
    end
end