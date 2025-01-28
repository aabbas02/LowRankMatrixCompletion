function krum_gradient = krum(stored_delta, L_byz)
    % Input:
    % stored_delta - a 3D array of size (n, r, L) containing the gradients
    % L - total number of gradients
    % L_byz - number of Byzantine gradients
    
    % Get dimensions
    [n, r, L] = size(stored_delta);
    
    % Reshape each gradient into a vector of size (n * r)
    gradients = reshape(stored_delta, n * r, L);
    
    % Initialize an array to store the scores for each vector
    scores = zeros(1, L);
    
    % Loop over each gradient to compute its score
    for i = 1:L
        % Compute squared Euclidean distances to all other vectors
        distances = zeros(1, L);
        for j = 1:L
            if i ~= j
                distances(j) = norm(gradients(:, i) - gradients(:, j))^2;
            else
                distances(j) = Inf; % set self-distance to Inf to ignore
            end
        end
        
        % Find the L - L_byz - 2 closest gradients to the i-th gradient
        [~, idx] = sort(distances); % sort distances
        closest_indices = idx(1:(L - L_byz - 2)); % get closest indices
        
        % Sum the distances to the closest vectors
        scores(i) = sum(distances(closest_indices));
    end
    
    % Find the gradient with the minimum score
    [~, min_index] = min(scores);
    
    % Reshape the selected gradient back to (n, r)
    krum_gradient = reshape(gradients(:, min_index), [n, r]);
end
