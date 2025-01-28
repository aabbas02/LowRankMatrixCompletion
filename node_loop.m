function Grad_U=node_loop(obj,U_t_prev,q,n,r)
            Grad_U=zeros(n,r);
            for k=1:q
               Bhat = ( obj.A_l(:, :, k) * U_t_prev ) \ obj.Y_l(:, k);
               Grad_U = Grad_U + obj.A_l(:, :, k)' * (obj.A_l(:, :, k) * U_t_prev * Bhat - obj.Y_l(:, k)) * Bhat';
            end
       end