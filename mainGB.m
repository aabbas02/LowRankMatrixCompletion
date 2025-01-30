% Parameters
tic
n = 200; % number of rows
q = 1000;  % number of columns
r = 4;  % rank
p = 0.4; % probability of observing each entry
Tgm = 10000;
mu = 2;  % parameter for projection operator
attck = 0;
L = 10;   % number of nodes
L_byz = 2;
C1 = 1;
% Generate Ustar, Bstar, and Xstar outside the function
Ustar = orth(randn(n, r));
Bstar = randn(r, q);
%G_B1 = 8;
G_B1 = 6;
G_B2 = 4;
G_B3 = 1;
qtilde = q/L;
% Randomly select 3 unique subset indices
% selected_indices = randperm(L, 3);
selected_indices = [1, 2, 3];
% Modify Bstar directly by multiplying selected subsets by G_B
Bstar1 = Bstar;
Bstar2 = Bstar;
Bstar3 = Bstar;
%Bstar4 = Bstar;
for idx = selected_indices
    % Define the column range for the current subset
    start_col = (idx-1) * qtilde + 1;
    end_col = idx * qtilde;
    % Multiply the subset by G_B
    Bstar1(:, start_col:end_col) = G_B1 * Bstar(:, start_col:end_col);
    Bstar2(:, start_col:end_col) = G_B2 * Bstar(:, start_col:end_col);
    Bstar3(:, start_col:end_col) = G_B3 * Bstar(:, start_col:end_col);
    %Bstar4(:, start_col:end_col) = G_B4 * Bstar(:, start_col:end_col);
end

Xstar1 = Ustar*Bstar1;
Xstar2 = Ustar*Bstar2;
Xstar3 = Ustar*Bstar3;

MC = 9;
error = zeros(3,Tgm,MC);
for mc = 1 : MC
    %----------------------------------
    %GB1
    %----------------------------------
    [~,sigma,~] = svds(Xstar1,r);
    kappa_tilde = sigma(1,1)/sigma(r,r);
    % Generate Y and S
    [Y,Ycol,rowIdx,Ycol_,rowsJ] = genY_and_IdxNew(Xstar1, p, L);
    [U0_init,S1,~] = svds(Y/p,r);
    U = U0_init(:,1:r);
    eta = 1/(S1(1,1)^2*p);    
    U_t_prev = U;
    nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    %stored_delta=zeros(n,r,L);
    for t = 1 : Tgm
        stored_delta = node_loopNewWithL(U_t_prev,q,n,r,rowIdx,Ycol,L,rowsJ,Ycol_);        
        byz_rev=byzantine_rev(L_byz,stored_delta,C1,attck);
        for j = 1 : L_byz
            idx=randi([1,L]);
            stored_delta(:,:,idx)=byz_rev(:,:,j);
        end
        I_t = [];
        for i = 1 : L
            U_temp = U_t_prev - eta * stored_delta(:,:,i);
            nrmsTemp = sqrt(sum(U_temp.^2,2)); %norm(U_temp')
            accept_update = all(nrmsTemp < (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 1.4 * mu * sqrt(r / n));
            % If all conditions are met, accept the update for this node
            if accept_update
                I_t = [I_t, i];  % Add index i to the set I_t
            else
                %fprintf('Node %d rejected based on norm condition.\n', i);
            end
         end
         accepted_deltas = stored_delta(:,:,I_t);
         Grad_U = krum(accepted_deltas, L_byz);
         U_cap=U_t_prev - eta * Grad_U;
         [Q1,~] = qr(U_cap);
         U_t_prev=Q1(:,(1:r));
         nrmsPrev = sqrt(sum(U_t_prev.^2,2));
         error(1,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
         if mod(t,100) == 0
             fprintf('GB1. Krum. Iteration %d, SD Error: %f\n', t, error(1,t,mc));
         end
    end
    %----------------------------------
    %GB2
    %----------------------------------
    [~,sigma,~] = svds(Xstar2,r);
    kappa_tilde = sigma(1,1)/sigma(r,r);
    % Generate Y and S
    [Y,Ycol,rowIdx,Ycol_,rowsJ] = genY_and_IdxNew(Xstar2, p, L);
    [U0_init,S1,~] = svds(Y/p,r);
    eta = 1/(S1(1,1)^2*p);    
    U = U0_init(:,1:r);
    U_t_prev = U;
    nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    %stored_delta=zeros(n,r,L);
    for t = 1 : Tgm
        stored_delta = node_loopNewWithL(U_t_prev,q,n,r,rowIdx,Ycol,L,rowsJ,Ycol_);        
        byz_rev=byzantine_rev(L_byz,stored_delta,C1,attck);
        for j = 1 : L_byz
            idx=randi([1,L]);
            stored_delta(:,:,idx)=byz_rev(:,:,j);
        end
        I_t = [];
        for i = 1 : L
            U_temp = U_t_prev - eta * stored_delta(:,:,i);
            nrmsTemp = sqrt(sum(U_temp.^2,2)); %norm(U_temp')
            accept_update = all(nrmsTemp < (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 1.4 * mu * sqrt(r / n));
            % If all conditions are met, accept the update for this node
            if accept_update
                I_t = [I_t, i];  % Add index i to the set I_t
            else
                %fprintf('Node %d rejected based on norm condition.\n', i);
            end
         end
         accepted_deltas = stored_delta(:,:,I_t);
         Grad_U = krum(accepted_deltas, L_byz);
         U_cap=U_t_prev - eta * Grad_U;
         [Q1,~] = qr(U_cap);
         U_t_prev=Q1(:,(1:r));
         nrmsPrev = sqrt(sum(U_t_prev.^2,2));
         error(2,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
         if mod(t,100) == 0
             fprintf('GB2. Krum. Iteration %d, SD Error: %f\n', t, error(2,t,mc));
         end
    end
    %----------------------------------
    %GB3
    %----------------------------------
    [~,sigma,~] = svds(Xstar3,r);
    kappa_tilde = sigma(1,1)/sigma(r,r);
    % Generate Y and S
    [Y,Ycol,rowIdx,Ycol_,rowsJ] = genY_and_IdxNew(Xstar3, p, L);
    [U0_init,S1,~] = svds(Y/p,r);
    eta = 1/(S1(1,1)^2*p);    
    U = U0_init(:,1:r);
    U_t_prev = U;
    nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    %stored_delta=zeros(n,r,L);
    for t = 1 : Tgm
        stored_delta = node_loopNewWithL(U_t_prev,q,n,r,rowIdx,Ycol,L,rowsJ,Ycol_);        
        byz_rev=byzantine_rev(L_byz,stored_delta,C1,attck);
        for j = 1 : L_byz
            idx=randi([1,L]);
            stored_delta(:,:,idx)=byz_rev(:,:,j);
        end
        I_t = [];
        for i = 1 : L
            U_temp = U_t_prev - eta * stored_delta(:,:,i);
            nrmsTemp = sqrt(sum(U_temp.^2,2)); %norm(U_temp')
            accept_update = all(nrmsTemp < (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 1.4 * mu * sqrt(r / n));
            % If all conditions are met, accept the update for this node
            if accept_update
                I_t = [I_t, i];  % Add index i to the set I_t
            else
                %fprintf('Node %d rejected based on norm condition.\n', i);
            end
         end
         accepted_deltas = stored_delta(:,:,I_t);
         Grad_U = krum(accepted_deltas, L_byz);
         U_cap=U_t_prev - eta * Grad_U;
         [Q1,~] = qr(U_cap);
         U_t_prev=Q1(:,(1:r));
         nrmsPrev = sqrt(sum(U_t_prev.^2,2));
         error(3,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
         if mod(t,100) == 0
             fprintf('GB3. Krum. Iteration %d, SD Error: %f\n', t, error(3,t,mc));
         end
    end

    mc
end

error_mean_MC = mean(error,3);
plotErrorMeanGB(error_mean_MC, n, q, r, p, L, L_byz, MC, C1, attck,G_B1,G_B2,G_B3)
toc
%------
