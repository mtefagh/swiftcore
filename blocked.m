function result = blocked(S, rev, solver)
% the currently available options for the LP solver are 'gurobi' and 'linprog'
    [m, n] = size(S);
    irev = m + find(rev == 0);
    model.obj = zeros(m+n, 1);
    model.obj(irev) = 1;
    model.A = [S.', -speye(n)];
    model.sense = repmat('=', n, 1);
    model.sense(rev == 0) = '<';
    model.rhs = zeros(n, 1);
    model.lb = [-Inf(m, 1); zeros(n, 1)];
    model.lb(irev) = -1;
    model.ub = [Inf(m, 1); zeros(n, 1)];
    if strcmp(solver, 'gurobi')
        params.outputflag = 0;
        result = gurobi(model, params);
        if ~strcmp(result.status, 'OPTIMAL')
            warning('Optimization is unstable!');
            fprintf('Optimization returned status: %s\n', result.status);
        end
    elseif strcmp(solver, 'linprog')
        problem.f = model.obj;
        problem.Aineq = model.A(rev == 0, :);
        problem.bineq = model.rhs(rev == 0);
        problem.Aeq = model.A(rev ~= 0, :);
        problem.beq = model.rhs(rev ~= 0);
        problem.lb = model.lb;
        problem.ub = model.ub;
        problem.solver = 'linprog';
        problem.options = optimset('Display', 'off');
        [result.x, result.objval, result.status, ~] = linprog(problem);
        if result.status ~= 1
            warning('Optimization is unstable!');
            fprintf('Optimization returned status: %s\n', result.status);
        end
    else
        model.b = model.rhs;
        model.c = model.obj;
        model.osense = 1;
        model.sense(model.sense == '=') = 'E';
        model.sense(model.sense == '<') = 'L';
        model.csense = model.sense;
        solution = solveCobraLP(model, 'solver', solver);
        result.x = solution.full;
        result.objval = solution.obj;
        result.status = solution.stat;
        if result.status ~= 1
            warning('Optimization is unstable!');
            fprintf('Optimization returned status: %s\n', result.status);
        end
    end
end