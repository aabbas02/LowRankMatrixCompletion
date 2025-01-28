% Parameters
n = 1000; % number of rows
q = 100;  % number of columns
r = 4;  % rank
L = 10;   % number of nodes
p = 0.1; % probability of observing each entry
T = 300;  % number of iterations
mu = 2;  % parameter for projection operator
L_byz = 25;
C1= 10;
m=30;
% Generate Ustar, Bstar, and Xstar outside the function
Ustar = orth(randn(n, r));
Bstar = randn(r, q);
Xstar = Ustar * Bstar;
[~,sigma,~] = svds(Xstar,r);
kappa_tilde = sigma(1,1)/sigma(r,r);
MC = 1;
error = zeros(4,T,MC);

for mc=1:MC
    % Generate Y and S 
    [Y, S] = generate_Y_and_S(Xstar, p);
    %[S, Y] = Generate(n,m,q,Xstar);

    [U0_init,S1,~] = svds(Y/p,r);
    U = U0_init(:,1:r);
    eta = 0.5/(S1(1,1)^2*p); 
    U_t_prev = U;
    
    %U = init(S,Y,r);
    %eta=0.9/(m*q);
    

    nodes=node.empty;
    for i=1:L
        n1=node;
        n1.id=i;
        n1.A_l=S(:,:,(i-1)*(q/L)+1:i*(q/L));
        n1.Y_l=Y(:,(i-1)*(q/L)+1:i*(q/L));
        nodes(i)=n1;
    end

    qtilde = q/L;
    
    % GD step Sum - Standard AltGdMin for LRMC - no robustness
    stored_delta=zeros(n,r,L);
    for t=1:T 
        for i=1:L
           stored_delta(:,:,i)=node_loop(nodes(i),U_t_prev,qtilde,n,r);
        end
        Grad_U = sum(stored_delta, 3);
        U_cap=U_t_prev - eta * Grad_U;
        [Q1,~] = qr(U_cap);
        U_t_prev=Q1(:,(1:r));
        error(1,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
        if mod(t,25) == 0
            fprintf('Sum. Iteration %d, SD Error: %f\n', t, error(1,t,mc));
        end
    end
    
    % GD step element wise median
    U_t_prev = U;
    %nrmsPrev
    stored_delta=zeros(n,r,L);
    for t=1:T 
        for i=1:L
           stored_delta(:,:,i)=node_loop(nodes(i),U_t_prev,qtilde,n,r);
        end
        byz_rev=byzantine_rev(L_byz,stored_delta,C1);
        %byz_rev=byzantine(L_byz,n,r);
        for j=1:L_byz
            idx=randi([1,L]);
            stored_delta(:,:,idx)=byz_rev(:,:,j);
        end
        I_t = [];
        for i=1:L
            U_temp = U_t_prev - eta * stored_delta(:,:,i);

           % Initialize a flag to check if all row conditions are met
           accept_update = true; % comment this

           % Check the norm condition for each row j in U_temp
           for j = 1:n
              % Calculate norms for condition check
              norm_U_temp_j = norm(U_temp(j, :));
              norm_U_t_prev_j = norm(U_t_prev(j, :));

              % Condition to check
              if norm_U_temp_j > (1 - 0.4 / (kappa_tilde^2)) * norm_U_t_prev_j + 1.4 * mu * sqrt(r / n)
                 accept_update = false;
                 break;
              end
           end
           %nrmTemp = norm(U_temp') % OR nrmTemp = sqrt(sum(U_temp.^2,2))
           %accept_update = all(nrmTemp > (1 - 0.4 / (kappa_tilde^2)) *
           %nrmPrev + 1.4 * mu * sqrt(r / n))
           %nrmPrev = nrmTemp

           % If all conditions are met, accept the update for this node
           if accept_update
              I_t = [I_t, i];  % Add index i to the set I_t
           else
              fprintf('Node %d rejected based on norm condition.\n', i);
           end   
        end
        accepted_deltas = stored_delta(:,:,I_t);
        Grad_U = coordinate_wise_median(accepted_deltas);
        U_cap=U_t_prev - eta * Grad_U;
        [Q1,~] = qr(U_cap);
        U_t_prev=Q1(:,(1:r));
        error(2,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
        if mod (t,25) == 0
            fprintf('Element. Iteration %d, SD Error: %f\n', t, error(2,t,mc));
        end
    end

    % GD step GM with check and attack
    U_t_prev = U;
    stored_delta=zeros(n,r,L);
    for t=1:T 
        for i=1:L
           stored_delta(:,:,i)=node_loop(nodes(i),U_t_prev,qtilde,n,r);
        end
        byz_rev=byzantine_rev(L_byz,stored_delta,C1);
        %byz_rev=byzantine(L_byz,n,r);
        for j=1:L_byz
            idx=randi([1,L]);
            stored_delta(:,:,idx)=byz_rev(:,:,j);
        end
        I_t = [];
        for i=1:L
            U_temp = U_t_prev - eta * stored_delta(:,:,i);

           % Initialize a flag to check if all row conditions are met
           accept_update = true;

           % Check the norm condition for each row j in U_temp
           for j = 1:n
              % Calculate norms for condition check
              norm_U_temp_j = norm(U_temp(j, :));
              norm_U_t_prev_j = norm(U_t_prev(j, :));

              % Condition to check
              if norm_U_temp_j > (1 - 0.4 / (kappa_tilde^2)) * norm_U_t_prev_j + 1.4 * mu * sqrt(r / n)
                 accept_update = false;
                 break;
              end
           end

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
        error(3,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
        if mod (t,25) == 0
             fprintf('GM. Iteration %d, SD Error: %f\n', t, error(3,t,mc));
        end
      
    end



    % GD step Krum with check and attack
    U_t_prev = U;
    stored_delta=zeros(n,r,L);
    for t=1:T 
        for i=1:L
           stored_delta(:,:,i)=node_loop(nodes(i),U_t_prev,qtilde,n,r);
        end
        byz_rev=byzantine_rev(L_byz,stored_delta,C1);
        %byz_rev=byzantine(L_byz,n,r);
        for j=1:L_byz
            idx=randi([1,L]);
            stored_delta(:,:,idx)=byz_rev(:,:,j);
        end
        I_t = [];
        for i=1:L
            U_temp = U_t_prev - eta * stored_delta(:,:,i);

           % Initialize a flag to check if all row conditions are met
           accept_update = true;

           % Check the norm condition for each row j in U_temp
           for j = 1:n
              % Calculate norms for condition check
              norm_U_temp_j = norm(U_temp(j, :));
              norm_U_t_prev_j = norm(U_t_prev(j, :));

              % Condition to check
              if norm_U_temp_j > (1 - 0.4 / (kappa_tilde^2)) * norm_U_t_prev_j + 1.4 * mu * sqrt(r / n)
                 accept_update = false;
                 break;
              end
           end

           % If all conditions are met, accept the update for this node
           if accept_update
              I_t = [I_t, i];  % Add index i to the set I_t
           else
              fprintf('Node %d rejected based on norm condition.\n', i);
           end   
        end
        accepted_deltas = stored_delta(:,:,I_t);
        Grad_U = krum(accepted_deltas, L_byz);
        U_cap=U_t_prev - eta * Grad_U;
        [Q1,~] = qr(U_cap);
        U_t_prev=Q1(:,(1:r));
        error(4,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
        if mod(t,25)==0
            fprintf('Krum. Iteration %d, SD Error: %f\n', t, error(4,t,mc));
        end
    end

    
    
end




% GD step GM
% stored_delta=zeros(n,r,L);
% for t=1:T 
%     for i=1:L
%        stored_delta(:,:,i)=node_loop(nodes(i),U_t_prev,qtilde,n,r);
%     end
%     Grad_U = do_gm_vec(stored_delta);
%     U_cap=U_t_prev - eta * Grad_U;
%     [Q1,~] = qr(U_cap);
%     U_t_prev=Q1(:,(1:r));
%     error(1,t) = subspace(Ustar,U_t_prev,n)/sqrt(r);
%     fprintf('Iteration %d, SD Error: %f\n', t, error(1,t));
% end

% GD step GM with check w/o attack
% U_t_prev = U;
% stored_delta=zeros(n,r,L);
% for t=1:T 
%     for i=1:L
%        stored_delta(:,:,i)=node_loop(nodes(i),U_t_prev,qtilde,n,r);
%     end
% %     byz_rev=byzantine_rev(L_byz,stored_delta,C1);
% %     for j=1:L_byz
% %         idx=randi([1,L]);
% %         stored_delta(:,:,idx)=byz_rev(:,:,j);
% %     end
%     I_t = [];
%     for i=1:L
%         U_temp = U_t_prev - eta * stored_delta(:,:,i);
%     
%        % Initialize a flag to check if all row conditions are met
%        accept_update = true;
%     
%        % Check the norm condition for each row j in U_temp
%        for j = 1:n
%           % Calculate norms for condition check
%           norm_U_temp_j = norm(U_temp(j, :));
%           norm_U_t_prev_j = norm(U_t_prev(j, :));
%         
%           % Condition to check
%           if norm_U_temp_j > (1 - 0.4 / (kappa_tilde^2)) * norm_U_t_prev_j + 1.4 * mu * sqrt(r / n)
%              accept_update = false;
%              break;
%           end
%        end
%     
%        % If all conditions are met, accept the update for this node
%        if accept_update
%           I_t = [I_t, i];  % Add index i to the set I_t
%        else
%           fprintf('Node %d rejected based on norm condition.\n', i);
%        end   
%     end
%     accepted_deltas = stored_delta(:,:,I_t);
%     Grad_U = do_gm_vec(accepted_deltas);
%     U_cap=U_t_prev - eta * Grad_U;
%     [Q1,~] = qr(U_cap);
%     U_t_prev=Q1(:,(1:r));
%     error(2,t) = subspace(Ustar,U_t_prev,n)/sqrt(r);
%     fprintf('Iteration %d, SD Error: %f\n', t, error(1,t));
% end

% Basic AltGDmin
% for t=1:T 
%     Grad_U=zeros(n,r);
%             for k=1:q
%                Bhat = ( S(:, :, k) * U_t_prev ) \ Y(:, k);
%                Grad_U = Grad_U + S(:, :, k)' * (S(:, :, k) * U_t_prev * Bhat - Y(:, k)) * Bhat';
%             end
% U_cap=U_t_prev - eta * Grad_U;
% [Q1,~] = qr(U_cap);
% U_t_prev=Q1(:,(1:r));
% error(1,t) = subspace(Ustar,U_t_prev,n)/sqrt(r);
% fprintf('Iteration %d, SD Error: %f\n', t, error(1,t));
% end
%%
% Plot SD error versus iteration
error_mean_MC = mean(error,3);
T1=1:T;
figure;
semilogy(T1,error_mean_MC(1,:), '->','MarkerSize',10, 'Color', 'black', 'MarkerIndices', 1:25:length(T1));
xlabel('Iteration')
ylabel('Error')
title('SD Error vs Iteration')
hold on 
semilogy(T1,error_mean_MC(2,:), '-diamond','MarkerSize',10, 'Color', 'blue', 'MarkerIndices', 1:25:length(T1));
semilogy(T1,error_mean_MC(3,:), '-*','MarkerSize',10, 'Color', '#A2142F', 'MarkerIndices', 1:25:length(T1));
semilogy(T1,error_mean_MC(4,:), '-o','MarkerSize',10, 'Color', 'red', 'MarkerIndices', 1:25:length(T1));
legend('Sum (No Attack)','Element-wise Median','GM','Krum')
hold off


% figure;
% semilogy(1:T, error, '-o');
% xlabel('Iteration');
% ylabel('SD Error');
% title('SD Error vs Iteration');
% grid on;


