% Parameters
tic
n = 1000; % number of rows
q = 500;  % number of columns
r = 3;  % rank
p = 0.2; % probability of observing each entry
Telem = 5000;  % number of iterations
Tgm = 5000;
mu = 2;  % parameter for projection operator
attck = 1;
L = 20;   % number of nodes
L_byz = 8;
C1 = 1;
% Generate Ustar, Bstar, and Xstar outside the function
Ustar = orth(randn(n, r));
Bstar = randn(r, q);
Xstar = Ustar * Bstar;
[~,sigma,~] = svds(Xstar,r);
kappa_tilde = sigma(1,1)/sigma(r,r);
MC = 101;
error = zeros(4,Telem,MC);
for mc = 1 : MC
    % Generate Y and S
    [Y,Ycol,rowIdx,Ycol_,rowsJ] = genY_and_IdxNew(Xstar, p, L);
    [U0_init,S1,~] = svds(Y/p,r);
    U = U0_init(:,1:r);
    eta = 1/(S1(1,1)^2*p);
    %--------------------------------------------------------
    % GD step element wise median
    %----------------------------------------------------------
    U_t_prev = U;
    nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    %stored_delta=zeros(n,r,L);
    for t = 1 : Telem
        %for i=1:L
           %stored_delta(:,:,i)=node_loop(nodes(i),U_t_prev,qtilde,n,r);
        %end
        stored_delta = node_loopNewWithL(U_t_prev,q,n,r,rowIdx,Ycol,L,rowsJ,Ycol_);
        byz_rev=byzantine_rev(L_byz,stored_delta,C1,attck);
        for j = 1 : L_byz
            idx=randi([1,L]);
            stored_delta(:,:,idx)=byz_rev(:,:,j);
        end
        I_t = [];
        for i = 1 : L
           U_temp = U_t_prev - eta * stored_delta(:,:,i);
           %{
           % Initialize a flag to check if all row conditions are met
           %accept_update = true; % comment this
           % Check the norm condition for each row j in U_temp
           %for j = 1:n
              % Calculate norms for condition check
              %norm_U_temp_j = norm(U_temp(j, :));
              %norm_U_t_prev_j = norm(U_t_prev(j, :));
              % Condition to check
              %if norm_U_temp_j > (1 - 0.4 / (kappa_tilde^2)) *
              %norm_U_t_prev_j + 1.4 * mu * sqrt(r / n) -
                 %accept_update = false;
                 %break;
              %end
           %end
           %}
           nrmsTemp = sqrt(sum(U_temp.^2,2)); %norm(U_temp')
           accept_update = all(nrmsTemp <= (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 1.4 * mu * sqrt(r / n));
           % If all conditions are met, accept the update for this node
           if accept_update
              I_t = [I_t, i];  % Add index i to the set I_t
           else
              %fprintf('Node %d rejected based on norm condition.\n', i);
           end   
        end
        accepted_deltas = stored_delta(:,:,I_t);
        Grad_U = coordinate_wise_median(accepted_deltas);
        U_cap = U_t_prev - eta * Grad_U;
        [Q1,~] = qr(U_cap);
        U_t_prev=Q1(:,(1:r));
        nrmsPrev = sqrt(sum(U_t_prev.^2,2));
        error(2,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
        if mod (t,100) == 0
            fprintf('Elmnt Median. n = %d, q = %d, r = %d, p = %f, L = %d, L_byz = %d, C1 = %d, Attack = %d. Iter. %d, SD Error: %f\n', n,q,r,p,L,L_byz,C1, attck, t, error(2,t,mc));
        end
    end
    %{
    %---------------------------------------------------------
    % GD step Sum - Standard AltGdMin for LRMC - no robustness
    %---------------------------------------------------------
    %stored_delta=zeros(n,r,L);
    %for t=1:T
    %    for i=1:L
    %       stored_delta(:,:,i)=node_loop(nodes(i),U_t_prev,qtilde,n,r);
    %    end
    %    Grad_U = sum(stored_delta, 3);
    %    U_cap=U_t_prev - eta * Grad_U;
    %    [Q1,~] = qr(U_cap);
    %    U_t_prev=Q1(:,(1:r));
    %    error(1,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
    %    if mod(t,25) == 0
    %        fprintf('Sum. Iteration %d, SD Error: %f\n', t, error(1,t,mc));
    %    end
    %end
    %}
    %---------------------------------
    % GD step GM with check and attacks
    %---------------------------------
    U_t_prev = U;
    nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    stored_delta=zeros(n,r,L);
    for t=1:1*Tgm
        %for i=1:L
           %stored_delta(:,:,i)=node_loopNew(U_t_prev,q,n,r,rowIdx,Xcol,rows,cols);
        %end
        stored_delta = node_loopNewWithL(U_t_prev,q,n,r,rowIdx,Ycol,L,rowsJ,Ycol_);        
        byz_rev=byzantine_rev(L_byz,stored_delta,C1,attck);
        for j=1:L_byz
            idx=randi([1,L]);
            stored_delta(:,:,idx)=byz_rev(:,:,j);
        end
        I_t = [];
        for i = 1 : L
            U_temp = U_t_prev - eta * stored_delta(:,:,i);
           %{
           % Initialize a flag to check if all row conditions are met
           %accept_update = true;
           % Check the norm condition for each row j in U_temp
           %for j = 1:n
           %   % Calculate norms for condition check
           %   norm_U_temp_j = norm(U_temp(j, :));
           %   norm_U_t_prev_j = norm(U_t_prev(j, :));
              % Condition to check
           %   if norm_U_temp_j > (1 - 0.4 / (kappa_tilde^2)) * norm_U_t_prev_j + 1.4 * mu * sqrt(r / n)
           %      accept_update = false;
           %      break;
           %   end
           %end
           %}
           nrmsTemp = sqrt(sum(U_temp.^2,2)); %norm(U_temp')
           accept_update = all(nrmsTemp <= (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 1.4 * mu * sqrt(r / n));
           % If all conditions are met, accept the update for this node
           if accept_update
              I_t = [I_t, i];  % Add index i to the set I_t
           else
              fprintf('Node %d rejected based on norm condition.\n', i);
           end   
        end
        accepted_deltas = stored_delta(:,:,I_t);
        Grad_U = do_gm_vec(accepted_deltas);
        U_cap=U_t_prev - eta * Grad_U;
        [Q1,~] = qr(U_cap);
        U_t_prev=Q1(:,(1:r));
        nrmsPrev = sqrt(sum(U_t_prev.^2,2));
        error(3,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
        if mod (t,100) == 0
             fprintf('GM.  n = %d, q = %d, r = %d, p = %f, L = %d, L_byz = %d, C1 = %d. Iteration %d, SD Error: %f\n', n,q,r,p, L,L_byz, C1, t, error(3,t,mc));
        end
    end
    %----------------------------------
    %GD step Krum with check and attack
    %----------------------------------
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
            %{
    %
    %        % Initialize a flag to check if all row conditions are met
    %        %accept_update = true;
    %
    %        % Check the norm condition for each row j in U_temp
    %        %for j = 1:n
    %           % Calculate norms for condition check
    %        %   norm_U_temp_j = norm(U_temp(j, :));
    %        %   norm_U_t_prev_j = norm(U_t_prev(j, :));
    %
    %           % Condition to check
    %        %   if norm_U_temp_j > (1 - 0.4 / (kappa_tilde^2)) * norm_U_t_prev_j + 1.4 * mu * sqrt(r / n)
    %        %      accept_update = false;
    %        %      break;
    %        %   end
    %        %end
            %}
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
         error(4,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
         if mod(t,100) == 0
             fprintf('Krum. Iteration %d, SD Error: %f\n', t, error(4,t,mc));
         end
    end
    mc
end
error_mean_MC = mean(error,3);
plotErrorMean(error_mean_MC, n, q, r, p, L, L_byz, MC, C1, attck)
toc
%------
