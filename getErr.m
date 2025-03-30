function [err] = getErr(U,r,q,numWrkrs,rowsJ,Xcol_,Y,idx)
    B = zeros(r,q);
    for j = 1 : numWrkrs
        tmp = zeros(r,q/numWrkrs);
        rowIdxj = rowsJ(j,:); 
        XcolJ = Xcol_(j,:);
        for k = 1 : q/numWrkrs
           rowIdx_jk = rowIdxj{k};
           tmp(:,k) = U(rowIdx_jk,:)\XcolJ{k};
        end
        offset = (j-1)*q/numWrkrs;        
        B(:,offset+1:(j*q)/numWrkrs) = tmp;
    end
    X = U*B;
    err = norm(X(idx) - Y(idx))/norm(Y(idx));
end
