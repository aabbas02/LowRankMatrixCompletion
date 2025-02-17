function [pi_map] = get_permutation_k(n,k)
    pi_map = 1:n;
    idx1 = randperm(n,k);
    %pi_map(idx1(randperm(length(idx1)))) = idx1;
    pi_map(idx1(randperm(k))) = idx1;
    pi_map = pi_map'; % transpose for column vector
end