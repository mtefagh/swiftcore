function reconstruction = swiftCore(S, rev, coreInd, varargin)
%% Usage
%  reconstruction = swiftCore(S, rev, coreInd [, solver])
%   * Inputs:
%   - S: the associated sparse stoichiometric matrix
%   - rev: the 0-1 vector with 1's corresponding to the reversible reactions
%   - coreInd: the set of indices corresponding to the core reactions
%   - solver: the LP solver to be used; the currently available options are
%   either 'gurobi' or 'linprog' with the default value of 'linprog'
%   * Outputs:
%   - reconstruction: the 0-1 indicator vector of the reactions constituting 
%   the consistent metabolic network reconstructed from the core reactions
    %% setting up the LP solver
    if ~isempty(varargin)
        solver = varargin{1};
    else
        solver = 'linprog';
    end
    [m, n] = size(S);
    reacNum = (1:n).';
    fullCouplings = (1:n).';
    %% finding the trivial full coupling relations
    flag = true;
    while flag
        flag = false;
        for i = m:-1:1
            if i <= size(S, 1)
                nzcols = find(S(i, :));
                % check to see if the i-th row of S has only two nonzero elements
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
                    fullCouplings(fullCouplings == reacNum(nzcols(2))) = reacNum(nzcols(1));
                    coreInd(coreInd == reacNum(nzcols(2))) = reacNum(nzcols(1));
                    reacNum(nzcols(2)) = [];
                    flag = true;
                end
            end
        end
    end
    S(:, rev == -1) = -S(:, rev == -1);
    rev(rev == -1) = 0;
    %% iterating the algorithm until the size of the subnetwork is no longer reduced by more than half
    originalCoreInd = ismember(reacNum, coreInd);
    coreInd = true(size(S, 2), 1);
    while true
        performance = sum(coreInd);
        coreInd = originalCoreInd;
        % the zero-tolerance parameter is the smallest flux value that is considered nonzero
        tol = norm(S(:, coreInd), 'fro')*eps(class(S));
        % identifying the blocked reversible reactions
        blocked = zeros(size(S, 2), 1);
        [Q, R, ~] = qr(transpose(S(:, coreInd)));
        Z = Q(:, sum(abs(diag(R)) > tol)+1:end);
        blocked(coreInd) = vecnorm(Z, 2, 2) < tol;
        % phase one of unblocking the reversible reactions
        c = 1;
        while any(blocked)
            % incrementing the core set until no reversible blocked reaction remains
            blockedSize = sum(blocked);
            [blocked, coreInd] = core(S, rev, blocked, coreInd, c, solver, tol);
            if 2*sum(blocked) > blockedSize
                c = 2*c;
            end
        end
        % phase two of unblocking the irreversible reactions
        [~, coreInd] = core(S, rev, blocked, coreInd, 0, solver, tol);
        S = S(:, coreInd);
        rev = rev(coreInd);
        originalCoreInd = originalCoreInd(coreInd);
        reacNum = reacNum(coreInd);
        if 2*sum(coreInd) > performance
            break;
        end
    end
    reconstruction = ismember(fullCouplings, reacNum);
end