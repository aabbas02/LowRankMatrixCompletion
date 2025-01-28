% Parameters
tic
n = 200; % number of rows
q = 1000;  % number of columns
r = 3;  % rank
L = 10;   % number of nodes
p = 0.3; % probability of observing each entry
Telem = 20000;  % number of iterations
Tgm = 20000;
mu = 2;  % parameter for projection operator
L_byz = 4;
C1= 1;
m=30;
% Generate Ustar, Bstar, and Xstar outside the function
Ustar = orth(randn(n, r));
Bstar = randn(r, q);
Xstar = Ustar * Bstar;
[~,sigma,~] = svds(Xstar,r);
kappa_tilde = sigma(1,1)/sigma(r,r);
MC = 1;
error = zeros(4,max(Telem,Tgm),MC);
for mc=1:MC
    % Generate Y and S
    [Y, S] = generate_Y_and_S(Xstar, p);
    %[S, Y] = Generate(n,m,q,Xstar);
    [U0_init,S1,~] = svds(Y/p,r);
    U = U0_init(:,1:r);
    eta = 1/(S1(1,1)^2*p);
    U_t_prev = U;
    
    nodes=node.empty;
    for i=1:L
        n1=node;
        n1.id=i;
        n1.A_l=S(:,:,(i-1)*(q/L)+1:i*(q/L));
        n1.Y_l=Y(:,(i-1)*(q/L)+1:i*(q/L));
        nodes(i)=n1;
    end
    qtilde = q/L;
    %--------------------------------------------------------
    % GD step element wise median
    %----------------------------------------------------------
    U_t_prev = U;
    nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    stored_delta=zeros(n,r,L);
    for t=1:1*Telem
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
           accept_update = all(nrmsTemp < (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 1.4 * mu * sqrt(r / n));
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
        nrmsPrev = sqrt(sum(U_t_prev.^2,2));
        error(2,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
        if mod (t,100) == 0
            fprintf('Element Wise Median. n = %d, q = %d, r = %d, p = %f, L = %d, L_byz = %d, C1 = %d. Iteration %d, SD Error: %f\n', n,q,r,p,L,L_byz,C1, t, error(2,t,mc));
        end
    end
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
    %---------------------------------
    % GD step GM with check and attacks
    %---------------------------------
    U_t_prev = U;
    nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    stored_delta=zeros(n,r,L);
    for t=1:1*Tgm
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
           accept_update = all(nrmsTemp < (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 1.4 * mu * sqrt(r / n));
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
    % % GD step Krum with check and attack
    % U_t_prev = U;
    % nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    % stored_delta=zeros(n,r,L);
    % for t=1:T
    %     for i=1:L
    %        stored_delta(:,:,i)=node_loop(nodes(i),U_t_prev,qtilde,n,r);
    %     end
    %     byz_rev=byzantine_rev(L_byz,stored_delta,C1);
    %     %byz_rev=byzantine(L_byz,n,r);
    %     for j=1:L_byz
    %         idx=randi([1,L]);
    %         stored_delta(:,:,idx)=byz_rev(:,:,j);
    %     end
    %     I_t = [];
    %     for i=1:L
    %         U_temp = U_t_prev - eta * stored_delta(:,:,i);
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
    %        nrmsTemp = sqrt(sum(U_temp.^2,2)); %norm(U_temp')
    %        accept_update = all(nrmsTemp < (1 - 0.4 / (kappa_tilde^2)) *nrmsPrev + 1.4 * mu * sqrt(r / n));
    %
    %        % If all conditions are met, accept the update for this node
    %        if accept_update
    %           I_t = [I_t, i];  % Add index i to the set I_t
    %        else
    %           fprintf('Node %d rejected based on norm condition.\n', i);
    %        end   
    %     end
    %     accepted_deltas = stored_delta(:,:,I_t);
    %     Grad_U = krum(accepted_deltas, L_byz);
    %     U_cap=U_t_prev - eta * Grad_U;
    %     [Q1,~] = qr(U_cap);
    %     U_t_prev=Q1(:,(1:r));
    %     nrmsPrev = sqrt(sum(U_t_prev.^2,2));
    %     error(4,t,mc) = subspace(Ustar,U_t_prev,n)/sqrt(r);
    %     if mod(t,25) == 0
    %         fprintf('Krum. Iteration %d, SD Error: %f\n', t, error(4,t,mc));
    %     end
    % end
    
    mc
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
figure;
semilogy(error_mean_MC(1,:), '->','MarkerSize',5, 'Color', 'black');
xlabel('Iteration')
ylabel('Error')
title("n = " + n + ", q = " + q +...
        ", r = " + r + ", p = " + p + ", L =" + L + ", $L_{byz} = $" + L_byz + '.',  ...
        'Interpreter', 'Latex', 'FontSize', 12');
hold on
semilogy(error_mean_MC(2,:), '-diamond','MarkerSize',5, 'Color', 'blue');
semilogy(error_mean_MC(3,:), '-*','MarkerSize',5, 'Color', '#A2142F');
semilogy(error_mean_MC(4,:), '-o','MarkerSize',5, 'Color', 'red');
legend('Sum (No Attack)','Element-wise Median','GM','Krum')
hold off
ID = (randi(10000));
stringTitle = ['MC_', num2str(MC), ...
               '_n_', num2str(n), '_q_', num2str(q), '_r_', ...
               num2str(r), '_L_', num2str(L), 'p', num2str(p),'_Lbyz_',num2str(L_byz),'C1_',num2str(C1),'_ID',num2str(ID)];
savefig([stringTitle, '.fig']);
save(['varsID_',num2str(ID),'.mat'], ...
    'n','q','r','p','error_mean_MC','MC','L','L_byz','C1','ID')
toc
% figure;
% semilogy(1:T, error, '-o');
% xlabel('Iteration');
% ylabel('SD Error');
% title('SD Error vs Iteration');
% grid on;
 