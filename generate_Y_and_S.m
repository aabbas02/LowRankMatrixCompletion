function [Y, S] = generate_Y_and_S(Xstar, p)
    % Generate the observed data matrix Y and selection matrices S
    % Inputs:
    %   - Xstar: original low-rank matrix
    %   - p: probability of observing each entry
    % Outputs:
    %   - Y: observed data matrix with unobserved entries set to zero
    %   - S: cell array containing diagonal matrices S_k for each column
    
    [n, q] = size(Xstar);
    Y = zeros(n, q);  % Initialize Y as a zero matrix
    S = zeros(n, n, q);   % Initialize cell array to store S_k matrices
    
    % Loop through each column to generate observed entries and S_k matrices
    for k = 1:q
        % Generate a binary vector indicating observed entries in column k
        observed_entries = rand(n, 1) < p;
        
        % Create the diagonal matrix S_k
        S(:,:,k) = diag(observed_entries);
        
        % Set observed entries in Y based on Xstar and S_k
        Y(:, k) = S(:,:,k) * Xstar(:, k);
    end
end
