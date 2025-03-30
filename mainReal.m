% Parameters
tic
clear all
close all
clc
r = 3;  % rank
Telem = 10000;  % number of iterations
Tgm = 10000;
mu = 2;  % parameter for projection operator
attck = 0;
L = 20;   % number of nodes
L_byz = 8;
C1 = 1e0;
MC = 1;
error = zeros(4,Telem,MC);
C_elem = 100; 
C_gm = 10;
C_krum = 100;
for mc = 1 : MC    
    [Ustar,Y,p] = getMovieLens(r);
    real = 1; idxFlag = 0; 
    [Y, rowIdx, ~, Ycol, ~,idx, Ycol_,rowsJ] = processMatrix(Y,p,real,idxFlag,L);
    %load('data10M.mat'); disp('done'); p = length(idx)/(n*q);
    n = size(Y,1); q = size(Y,2);    
    [U0_init,S1,~] = svds(Y/p,r);
    kappa_tilde = S1(1,1)/S1(r,r);    
    U = U0_init(:,1:r);
    eta_elem = C_elem/(S1(1,1)^2*p);
    eta_gm = C_gm/(S1(1,1)^2*p);
    eta_krum = C_krum/(S1(1,1)^2*p);
    
    %--------------------------------------------------------
    % GD step element wise median
    %----------------------------------------------------------
    U_t_prev = U;
    nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    for t = 1 : Telem
        [stored_delta] = nodeLoopReal(U_t_prev,q,n,r,rowIdx,Ycol,L,rowsJ,Ycol_);
        byz_rev = byzantine_rev(L_byz,stored_delta,C1,attck);
        for j = 1 : L_byz
            idx_ = randi([1,L]);
            stored_delta(:,:,idx_) = byz_rev(:,:,j);
        end
        I_t = [];
        for i = 1 : L
            U_temp = U_t_prev - eta_elem * stored_delta(:,:,i);
            nrmsTemp = sqrt(sum(U_temp.^2,2)); %norm(U_temp')
            accept_update = all(nrmsTemp <= (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 1.4 * mu * sqrt(r / n));
            if accept_update
                I_t = [I_t, i];  % Add index i to the set I_t
            else
                if mod(t,100) == 0
                    fprintf('Elem. Node %d rejected based on norm condition.\n', i);
                end
            end
        end
        accepted_deltas = stored_delta(:,:,I_t);
        Grad_U = coordinate_wise_median(accepted_deltas);
        U_cap = U_t_prev - eta_elem * Grad_U;
        [Q1,~] = qr(U_cap);
        U_t_prev=Q1(:,(1:r));
        nrmsPrev = sqrt(sum(U_t_prev.^2,2));
        error(2,t,mc) = getErr(U_t_prev,r,q,L,rowsJ,Ycol_,Y,idx);
        if mod (t,100) == 0
            fprintf('Elmnt Median. n = %d, q = %d, r = %d, p = %f, L = %d, L_byz = %d, C1 = %d, Attack = %d. Iter. %d, Rel. Error: %f\n', n,q,r,p,L,L_byz,C1, attck, t, error(2,t,mc));
        end
    end
    %---------------------------------
    % GD step GM with check and attacks
    %---------------------------------
    % U_t_prev = U;
    % nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    % for t = 1 : 1 * Tgm
    %     [stored_delta] = nodeLoopReal(U_t_prev,q,n,r,rowIdx,Ycol,L,rowsJ,Ycol_);
    %     byz_rev=byzantine_rev(L_byz,stored_delta,C1,attck);
    %     for j = 1 : L_byz
    %         idx_ = randi([1,L]);
    %         %idx = 1;
    %         stored_delta(:,:,idx_) = byz_rev(:,:,j);
    %     end
    %     I_t = [];
    %     for i = 1 : L
    %         U_temp = U_t_prev - eta_gm * stored_delta(:,:,i);
    %         nrmsTemp = sqrt(sum(U_temp.^2,2)); %norm(U_temp')
    %         accept_update = all(nrmsTemp <= (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 0.4 * mu * sqrt(r / n));
    %         % If all conditions are met, accept the update for this node
    %         if accept_update
    %             I_t = [I_t, i];  % Add index i to the set I_t
    %         else
    %             if mod(t,100) == 0
    %                 fprintf('GM. Node %d rejected based on norm condition.\n', i);
    %             end
    %         end
    %     end
    %     accepted_deltas = stored_delta(:,:,I_t);
    %     Grad_U = do_gm_vec(accepted_deltas);
    %     U_cap=U_t_prev - eta_gm * Grad_U;
    %     [Q1,~] = qr(U_cap);
    %     U_t_prev=Q1(:,(1:r));
    %     nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    %     error(3,t,mc) = getErr(U_t_prev,r,q,L,rowsJ,Ycol_,Y,idx);
    %     if mod (t,100) == 0
    %         fprintf('GM.  n = %d, q = %d, r = %d, p = %f, L = %d, L_byz = %d, C1 = %d. Iteration %d, Rel. Error: %f\n', n,q,r,p, L,L_byz, C1, t, error(3,t,mc));
    %     end
    % end
    %----------------------------------
    %GD step Krum with check and attack
    %----------------------------------
    U_t_prev = U;
    nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    for t = 1 : Tgm
        [stored_delta] = nodeLoopReal(U_t_prev,q,n,r,rowIdx,Ycol,L,rowsJ,Ycol_);
        byz_rev=byzantine_rev(L_byz,stored_delta,C1,attck);
        for j = 1 : L_byz
            idx_ = randi([1,L]);
            stored_delta(:,:,idx_) = byz_rev(:,:,j);
        end
        I_t = [];
        for i = 1 : L
            U_temp = U_t_prev - eta_krum * stored_delta(:,:,i);
            nrmsTemp = sqrt(sum(U_temp.^2,2)); %norm(U_temp')
            accept_update = all(nrmsTemp < (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 1.4 * mu * sqrt(r / n));
            % If all conditions are met, accept the update for this node
            if accept_update
                I_t = [I_t, i];  % Add index i to the set I_t
            else
                if mod(t,100) == 0
                    fprintf('Krum. Node %d rejected based on norm condition.\n', i);
                end
            end
        end
        accepted_deltas = stored_delta(:,:,I_t);
        Grad_U = krum(accepted_deltas, L_byz);
        U_cap=U_t_prev - eta_krum * Grad_U;
        [Q1,~] = qr(U_cap);
        U_t_prev=Q1(:,(1:r));
        nrmsPrev = sqrt(sum(U_t_prev.^2,2));
        error(4,t,mc) = getErr(U_t_prev,r,q,L,rowsJ,Ycol_,Y,idx);
        if mod(t,100) == 0
            fprintf('Krum. Iteration %d, Rel Error: %f\n', t, error(4,t,mc));
        end
    end
    mc
end
error_mean_MC = mean(error,3);
plotErrorMean(error_mean_MC, n, q, r, p, L, L_byz, MC, C1, attck,real,C_elem,C_gm,C_krum)
toc
