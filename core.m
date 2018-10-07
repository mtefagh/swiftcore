function flux = core(S, rev, blocked, weights, solver)
%% the currently available options for the LP solver are 'gurobi' and 'linprog'
    [m, n] = size(S);
    dense = zeros(n, 1);
    dense(blocked == 1) = normrnd(0, 1, [sum(blocked), 1]);
    coreInd = weights == 0;
    k = n - sum(coreInd);
    model.obj = [dense; weights(coreInd == 0)];
    temp = speye(n);
    model.A = [S,sparse(m,k); temp(coreInd==0,:),speye(k); -temp(coreInd==0,:),speye(k)];
    model.sense = repmat('=', m+2*k, 1);
    model.sense(m+1:m+2*k) = '>';
    model.rhs = zeros(m+2*k, 1);
    model.lb = -Inf(n, 1);
    model.lb(blocked == 1) = -1;
    model.lb(rev == 0) = 0;
    if ~any(blocked)
        model.lb(coreInd ~= 0 & rev == 0) = 1;
    end
    model.lb = [model.lb; -Inf(k, 1)];
    model.ub = Inf(n, 1);
    model.ub(blocked == 1) = 1;
    model.ub = [model.ub; Inf(k, 1)];
    if strcmp(solver, 'gurobi')
        params.outputflag = 0;
        result = gurobi(model, params);
        if strcmp(result.status, 'OPTIMAL')
            flux = result.x(1:n);
        else
            warning('Optimization was stopped with status %s\n', result.status);
        end
    elseif strcmp(solver, 'linprog')
        problem.f = model.obj;
        problem.Aineq = -model.A(m+1:m+2*k, :);
        problem.bineq = model.rhs(m+1:m+2*k);
        problem.Aeq = model.A(1:m, :);
        problem.beq = model.rhs(1:m);
        problem.lb = model.lb;
        problem.ub = model.ub;
        problem.solver = 'linprog';
        problem.options = optimset('Display', 'off');
        [result.x, result.objval, result.status, ~] = linprog(problem);
        if result.status == 1
            flux = result.x(1:n);
        else
            warning('Optimization was stopped with status %s\n', result.status);
        end
    else
        model.b = model.rhs;
        model.c = model.obj;
        model.osense = 1;
        model.sense(model.sense == '=') = 'E';
        model.sense(model.sense == '>') = 'G';
        model.csense = model.sense;
        solution = solveCobraLP(model, 'solver', solver);
        result.x = solution.full;
        result.objval = solution.obj;
        result.status = solution.stat;
        if result.status == 1
            flux = result.x(1:n);
        else
            warning('Optimization was stopped with status %s\n', result.status);
        end
    end
end