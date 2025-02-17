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