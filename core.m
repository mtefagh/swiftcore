function [blocked, coreInd] = core(S, rev, blocked, coreInd, c, solver, epsilon)
    [m, n] = size(S);
    denseSet = zeros(n, 1);
    denseSet(blocked == 1) = normrnd(0, c, [sum(blocked), 1]);
    sparseSet = ones(n, 1);
    sparseSet(coreInd) = 0;
    k = sum(sparseSet);
    model.obj = [denseSet; ones(k, 1)];
    temp = speye(n);
    model.A = [S, sparse(m,k); temp(sparseSet == 1, :), speye(k); -temp(sparseSet == 1, :), speye(k)];
    model.sense = repmat('=', m+2*k, 1);
    model.sense(m+1:m+2*k) = '>';
    model.rhs = zeros(m+2*k, 1);
    model.lb = -Inf(n, 1);
    model.lb(blocked == 1) = -1;
    model.lb(rev == 0) = 0;
    if  c == 0
        model.lb(sparseSet + rev == 0) = 1;
    end
    model.lb = [model.lb; -Inf(k, 1)];
    model.ub = Inf(n, 1);
    model.ub(blocked == 1) = 1;
    model.ub = [model.ub; Inf(k, 1)];
    if strcmp(solver, 'gurobi')
        params.outputflag = 0;
        result = gurobi(model, params);
        if strcmp(result.status, 'OPTIMAL')
            index = abs(result.x(1:n)) > epsilon;
            blocked(index) = 0;
            coreInd = max(coreInd, index);
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
            index = abs(result.x(1:n)) > epsilon;
            blocked(index) = 0;
            coreInd = max(coreInd, index);
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
            index = abs(result.x(1:n)) > epsilon;
            blocked(index) = 0;
            coreInd = max(coreInd, index);
        else
            warning('Optimization was stopped with status %s\n', result.status);
        end
    end
end