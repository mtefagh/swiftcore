function consistent = swiftcc(S, rev, varargin)
%% Usage
%  consistent = swiftcc(S, rev [, solver])
%   * Inputs:
%   - S: the associated sparse stoichiometric matrix
%   - rev: the 0-1 vector with 1's corresponding to the reversible reactions 
%   - solver: the LP solver to be used; the currently available options are
%   'gurobi', 'linprog', and 'cplex' with the default value of 'linprog'
%   * Outputs:
%   - consistent: the 0-1 indicator vector of the reactions constituting 
%   the maximum flux consistent metabolic subnetwork
    %% setting up the LP solver
    if ~isempty(varargin)
        solver = varargin{1};
    else
        solver = 'linprog';
    end
    [m, n] = size(S);
    consistent = true(n, 1);
    %% identifying the blocked irreversible reactions
    result = blocked(S, rev, solver);
    consistent(result.x(m+1:end) < -0.5) = false;
    %% setting up the zero-tolerance parameter
    tol = norm(S(:, consistent), 'fro')*eps(class(S));
    %% identifying the blocked reversible reactions
    [Q, R, ~] = qr(transpose(S(:, consistent)));
    Z = Q(:, sum(abs(diag(R)) > tol)+1:end);
    %% finding the consistent reactions of the original metabolic network
    consistent(consistent & rev == 1) = vecnorm(Z(rev(consistent) == 1, :), ...
        2, 2) > tol;
    consistent = find(consistent);
end