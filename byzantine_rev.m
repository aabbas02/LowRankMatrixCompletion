function byz_rev=byzantine_rev(L_byz,stored_delta,C1,attck)
   [size1, size2, ~] = size(stored_delta);
   % s = mean(stored_delta,3);
   s = sum(stored_delta,3);
   byz_rev=zeros(size1,size2,L_byz);
   for i = 1 : L_byz
      if attck 
        byz_rev(:,:,i) = -C1*s;
      else
        byz_rev(:,:,i) = 1 + randn(size1,size2);
      end
   end
   
end