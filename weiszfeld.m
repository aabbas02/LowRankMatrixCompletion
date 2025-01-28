
function [guess, iters]= weiszfeld(points)
  tol=1e-7;
  
    function distances=distance_func(x)
        distances=pdist2(x,points);
    end
   
% initial guess
% M = zeros(1,size(points,1));
% for i=1:size(points,1)
%     tmp = distance_func(points(i,:));
%     M(1,i)=sum(tmp,2);
% end
% 
% [~,I] = min(M);
% distances=distance_func(points(I,:));
% distances(distances==0)=1;
% S=sum((1./distances))-1;
% S1=zeros(1,size(points,2));
% for i=1:size(distances,2)
%     if i==I
%         break
%     end
%     S1=S1+(points(I,:)-points(i,:))./distances(1,i);
% end
% d_p = S1/norm(S1);
% t_p = (norm(S1)-1)/S;
% disp(S)
% guess = points(I,:) + t_p*d_p; 

guess=mean(points,1);
iters=0;
guess_movement=1;
while iters < 6000 && guess_movement >= tol
    distances=distance_func(guess);
    distances(distances==0)=1;
    S=sum((1./distances));
    S1=zeros(1,size(points,2));
    for i=1:size(distances,2)
        S1=S1+points(i,:)./distances(1,i);
    end
    guess_next = S1/S;
    guess_movement = norm(guess - guess_next);
    guess = guess_next;
    iters=iters+1;

end
end
