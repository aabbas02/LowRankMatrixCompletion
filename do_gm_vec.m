function agg1=do_gm_vec(node_compute)  
  [size1, size2 , L] = size(node_compute);
  aggr=zeros(L,size1*size2);
  for i=1:L
      aggr(i,:)=reshape(node_compute(:,:,i),[1,size1*size2]);
  end
  agg=weiszfeld(aggr);
  % Find the closest gradient to the aggregate 'agg' in terms of squared l2 distance
  distances = zeros(1, L);
  for i = 1:L
      distances(i) = norm(aggr(i, :) - agg)^2;
  end
  % Find the index of the closest gradient
  [~, min_index] = min(distances);
    
  % Reshape the selected closest gradient back to (size1, size2)
  agg1 = reshape(aggr(min_index, :), [size1, size2]);
  % agg1=reshape(agg,[size1,size2]);
end