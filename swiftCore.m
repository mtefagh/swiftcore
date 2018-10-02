function coreInd = swiftCore(S, rev, coreInd, solver)
    [m, n] = size(S);
    reacNum = (1:n).';
    fullCouplings = (1:n).';
    %% identifying the fully coupled pairs of reactions
    % finding the trivial full coupling relations
    flag = true;
    while flag
        flag = false;
        for i = m:-1:1
            if i <= size(S, 1)
                nzcols = find(S(i, :));
                % check to see if the i-th row of S has only 2 nonzero elements
                if length(nzcols) == 2
                    % deleting the reaction from the rev vector
                    if rev(nzcols(2)) ~= 1
                        if S(i, nzcols(1))/S(i, nzcols(2)) < 0
                            rev(nzcols(1)) = rev(nzcols(2));
                        else
                            rev(nzcols(1)) = -1 - rev(nzcols(2));
                        end
                    end
                    rev(nzcols(2)) = [];
                    % merging the fully coupled pair of reactions
                    S(:, nzcols(1)) = S(:, nzcols(1)) - S(i, nzcols(1))/S(i, nzcols(2))*S(:, nzcols(2));
                    S(:, nzcols(2)) = [];
                    % deleting the zero rows from the stoichiometric matrix
                    S = S(any(S, 2), :);
                    fullCouplings(fullCouplings == reacNum(nzcols(2))) = ...
                        reacNum(nzcols(1));
                    coreInd(coreInd == reacNum(nzcols(2))) = reacNum(nzcols(1));
                    reacNum(nzcols(2)) = [];
                    flag = true;
                end
            end
        end
    end
    S(:, rev == -1) = -S(:, rev == -1);
    rev(rev == -1) = 0;
    
    n = size(S, 2); 
    originalCoreInd = ismember(reacNum, coreInd);
    coreInd = true(n, 1);
    while true
        temp1 = sum(coreInd);
        coreInd = originalCoreInd;
        n = size(S, 2);
        blocked = zeros(n, 1);
        % identifying the blocked reversible reactions
        [Q, R, ~] = qr(transpose(S(:, coreInd)));
        % setting up the zero-tolerance parameter i.e. the smallest flux value 
        % that is considered nonzero
        tol = norm(S(:, coreInd), 'fro')*eps(class(S));
        Z = Q(:, sum(abs(diag(R)) > tol)+1:end);
        blocked(coreInd) = vecnorm(Z, 2, 2) < tol;
        c = 1;
        while any(blocked)
            temp2 = sum(blocked);
            [blocked, coreInd] = core(S, rev, blocked, coreInd, c, solver, tol);
            if 2*sum(blocked) > temp2
                c = 2*c;
            end
        end
        [~, coreInd] = core(S, rev, blocked, coreInd, 0, solver, tol);
        S = S(:, coreInd);
        rev = rev(coreInd);
        originalCoreInd = originalCoreInd(coreInd);
        reacNum = reacNum(coreInd);
        if 2*sum(coreInd) > temp1
            break;
        end
    end
    coreInd = ismember(fullCouplings, reacNum);
end