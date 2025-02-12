function [SDVals] = altGDMin_MtrxSensingPerm(Ak_, ykPerm_,AkCllps_,ykCllps_,Uinit,r,T,Ustr,r_,updtP,same)
	% This function should implement altGDMin with both permuted and non-permuted measurements
	% The above functionality is achieved by setting the arguments correspondingly
	% AltGDMin wout Perm: updtP = 0, Uinit = U0, Ak_ = Ak, ykPerm_ = yk, AkCllps_ = Ak, ykCllps_ = yk
	% AltGDMin with Perm: updtP = 1, Uinit = U0Cllps, Ak_ = Ak, ykPerm_ = ykPerm, AkCllps_ = AkCllps_, ykCllps_ = ykCllps_
	%---
    % Algorithm
    % Init: U^(0), B^(0)
    % Steps: min P, min U, min B
    % if t == 0, update bk by collapsed estimate, else update by full
    % estimate
    % repeat
    m = size(Ak_{1}, 1);
    n = size(Ak_{1}, 2);
    SDVals = zeros(T+1,1);
    U = Uinit;
    SDVals(1) = norm((eye(n) - U*U')*Ustr ,'fro');
    q = length(ykPerm_);
    B = zeros(r,q);
    gradU = zeros(n,r);
    %if updtP && same
        yHat = zeros(m,q);
        yPerm = zeros(m,q);
        for k =  1 : q
            yPerm(:,k) = ykPerm_{k};
        end
    %end
    for i = 1 : T
        % b_k update and P_k update in for k = 1 : q loop below
        for k = 1 : q
            % Least-squares B_k update
            if i == 1 % collapsed least-squares
                B(:,k) = (AkCllps_{k}*U)\ykCllps_{k};
            else % full measurements least-squares
                B(:,k) = (Ak_{k}*U)\ykPerm_{k};
            end
            if updtP 
                yHat(:,k) = Ak_{k}*U*B(:,k);
            end
            % min over P_k
            if updtP && same == 0
                for s = 1 : length(r_)
                    start = sum(r_(1:s)) - r_(s) + 1;
                    stop = sum(r_(1:s));
                    C = yPerm(start:stop,k)*yHat(start:stop,k)';                    
                    M = matchpairs(-C,1e10); % M is a matrix with 2 columns and m rows, 
                                             % The second column has ascending
                                             % indices in order 1, ..., m
                                             % The first column has the
                                             % corresponding/matching row indices
                                             % 5,1 means P(5,1) = 1, i.e., 
                                             % row 5 matched to 1
                    idx  = M(:,1);
                    idx = start - 1  + idx;
                    Ak_{k}(idx,:) = Ak_{k}(start:stop,:);                   
                    %[~,idx1] = sort(yHat_k(start:stop));
                    %[~,idx2] = sort(ykPerm_{k}(start:stop));
                    %idx1 = start - 1 + idx1;
                    %idx2 = start - 1 + idx2;
                    %Ak_{k}(idx2,:) = Ak_{k}(idx1,:);
                end
            end
        end % exit the for loop over columns 1 through q
        if updtP && same == 1
            for s =  1 : length(r_)
                start = sum(r_(1:s)) - r_(s) + 1;
                stop = sum(r_(1:s));
                C = yPerm(start:stop,:)*yHat(start:stop,:)';
                M = matchpairs(-C,1e100); % M is a matrix with 2 columns and m rows, 
                                         % The second column has ascending
                                         % indices in order 1, ..., m
                                         % The first column has the
                                         % corresponding/matching row indices
                                         % 5,1 means P(5,1) = 1, i.e., 
                                         % row 5 matched to 1
                idx  = M(:,1);
                idx = start - 1  + idx;
                for k = 1:q
                    Ak_{k}(idx,:) = Ak_{k}(start:stop,:);                   
                end
                %M(M(:,1)) = M(:,2);
                %temp = M(:,1);
                %temp = start - 1 + temp;
                %assignment(start:stop) = temp;   
            end
        end
        % U update
        X = U*B;
        if i == 1
            X0 = X;
        end
        gradU = 0*gradU;
        for k = 1 : q
            gradU = gradU + Ak_{k}'*(Ak_{k}*X(:,k)-ykPerm_{k})*B(:,k)';
        end
        eta = 5e-1/norm(X0)^2;
        U = U - (eta/m)*gradU;
        [U,~,~] = qr(U,'econ');
        SDVals(i + 1) = norm( (eye(n) - U*U')*Ustr ,'fro' );
    end
end 