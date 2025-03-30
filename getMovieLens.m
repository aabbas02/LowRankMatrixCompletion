%-------------
function [Ustr,X,p] = getMovieLens(r)
    A = readmatrix("ratings100K.xlsx");
    %--------------------
    %load("ratings1M.mat");
    %A = ratings1M;
    %--------------------
    %load("ratings10M.mat")
    n = max(A(:,1));
    q = max(A(:,2));
    X = zeros(n,q);
    num = 0;
    for k = 1 : size(A,1)
        i = A(k,1);
        j = A(k,2);
        X(i,j) = A(k,3);
        num = num + 1;
    end
    p = num/(n*q);
    %X = X';
    [Ustr,~,~] = svds(X,r);
    Ustr = Ustr(:,1:r);
end
%-------------
