function median_gradient = coordinate_wise_median(stored_delta)
    % Input:
    % stored_delta - a 3D array of size (n, r, L) containing the gradients
    
    % Get dimensions
    [n, r, L] = size(stored_delta);
    
    % Reshape each gradient into a vector of size (n * r)
    vectorized_gradients = reshape(stored_delta, n * r, L);
    
    % Calculate the coordinate-wise median across the L gradients
    coordinate_median = median(vectorized_gradients, 2);
    
    % Reshape the median result back to (n, r)
    median_gradient = reshape(coordinate_median, [n, r]);
end
