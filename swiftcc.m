function consistent = swiftcc(S, rev, varargin)
%% swiftcc finds the consistent reactions of the metabolic network
    % setting up the LP solver
    if ~isempty(varargin)
        solver = varargin{1};
    else
        solver = 'linprog';
    end
    [m, n] = size(S);
    consistent = true(n, 1);
    % identifying the blocked irreversible reactions
    result = blocked(S, rev, solver);
    consistent(result.x(m+1:end) < -0.5) = false;
    % setting up the zero-tolerance parameter
    tol = norm(S(:, consistent), 'fro')*eps(class(S));
    % identifying the blocked reversible reactions
    [Q, R, ~] = qr(transpose(S(:, consistent)));
    Z = Q(:, sum(abs(diag(R)) > tol)+1:end);
    consistent(consistent & rev == 1) = vecnorm(Z(rev(consistent) == 1, :), ...
        2, 2) > tol;
    consistent = find(consistent);
end