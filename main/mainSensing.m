clc
close all
clear all
dir = pwd;
% For linux, replace '\' with '/'
cd ..
addpath(genpath('.\functions'));
addpath(genpath('.\functionsMtrxSnsng'));
cd(dir)    
%---------------------------------
r = 5;
n = 600;
q = 1000;
m = 100;
numBlocks = 8;   %effectively, m_new = numBlocks
r_ = ones(1,numBlocks)*(m/numBlocks);
MC = 1;
% generate rank-r X*
Ustr = orth(randn(n,r));
Bstr = randn(r,q);
X  = Ustr*Bstr;
% generate q matrices A_k of size m x n, m << n
Ak_ = cell(q,1);
AkCllps_ = cell(q,1);
yk_ = cell(q,1);
ykCllps_ = cell(q,1);
ykPerm_ = cell(q,1);
M = zeros(n,q);
MCllps = zeros(n,q);
MPerm = zeros(n,q);
for k = 1 : q
    Ak_{k} = randn(m,n);
    yk_{k} = Ak_{k}*X(:,k);
    pi_map = get_permutation_r(m,r_);
    ykPerm_{k} = yk_{k}(pi_map);
    AkCllps_{k} = zeros(length(r_),n);
    ykCllps_{k} = zeros(length(r_),1);
        for i = 1 : length(r_)
            start = sum(r_(1:i)) - r_(i) + 1;
            stop  = sum(r_(1:i));
            AkCllps_{k}(i,:) = sum(Ak_{k}(start:stop,:));
            ykCllps_{k}(i) = sum(ykPerm_{k}(start:stop,:));
        end
    M(:, k)  = Ak_{k}'*yk_{k}; 
    MCllps(:, k) = AkCllps_{k}'*ykCllps_{k};
    MPerm(:, k) = Ak_{k}'*ykPerm_{k};
end
[U0,~,~,] = svd(M,"econ");
U0 = U0(:,1:r);
%-----------------------------------
[U0Cllps,~,~] = svd(MCllps,"econ");
U0Cllps = U0Cllps(:,1:r);
%---------------------------------------
[U0Perm,~,~] = svd(MPerm,"econ");
U0Perm = U0Perm(:,1:r);
%------------------------------
SDU0 = norm((eye(n) - Ustr*Ustr')*U0);
SDU0Cllps = norm((eye(n) -  Ustr*Ustr')*U0Cllps);
SDU0Perm = norm((eye(n) - Ustr*Ustr')*U0Perm);
T = 100;
SDVals_UnPerm = altGDMin_MtrxSensing(Ak_, yk_, U0,r,T,Ustr);
SDVals_sLcl = altGDMin_MtrxSensingPerm(Ak_, ykPerm_,AkCllps_, ykCllps_, U0Cllps,r,T,Ustr,r_);
SDVals_Perm = altGDMin_MtrxSensing(Ak_, ykPerm_, U0Perm,r,T,Ustr);

figure;
semilogy(SDVals_sLcl, ...
    'DisplayName', 'SDVals (s-Local)', 'LineWidth', 1.45, 'Marker', 'o', 'MarkerSize', 7);
hold on;
semilogy(SDVals_Perm, ...
    'DisplayName', 'SDVals (Naive)', 'LineWidth', 1.45, 'Marker', 'square', 'MarkerSize', 7);
%
semilogy(SDVals_UnPerm, ...
    'DisplayName', 'SDVals (Unpermuted)', 'LineWidth', 1.45, 'Marker', 'x', 'MarkerSize', 7);
grid on

stringTitle = ['MtrxSnsngPerm_MC_', num2str(MC), ...
               '_n_', num2str(n), '_q_', num2str(q), '_r_', num2str(r), '_m_', num2str(m)];

savefig([stringTitle, '.fig']);

title("n = " + n + ", q = " + q +...
      ", r = " + r + ", m = " + m + ", num Blocks = " + numBlocks,...
       'Interpreter', 'Latex', 'FontSize',14)

legend('Interpreter', 'Latex', 'Fontsize', 9);
ylabel("$SD(U,U^*)$","FontSize",14,'Interpreter','Latex')
xlabel('Iterations (t)', 'FontSize',14, 'Interpreter','Latex')
function [pi_map] = get_permutation_r(n,r_)
    pi_map = zeros(n,1);
    for t = 1 : length(r_)
        start  = sum(r_(1:t)) - r_(t) + 1;
        stop   = sum(r_(1:t));
        idx    = start:stop;
        idx    = idx(randperm(length(idx)));
        pi_map(start:stop) = idx;
    end
end
