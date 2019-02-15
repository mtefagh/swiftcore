function component = partition(model, solver, algorithm)
    S = model.S;
    rev = model.rev;
    lb = model.lb;
    ub = model.ub;
    c = model.c;
    rxns = model.rxns;
    [m, n] = size(S);
    DG = zeros(m+1);
    for i = 1:n
        head = S(:, i) > 0;
        if ~any(head)
            head = m+1;
        end
        tail = S(:, i) < 0;
        if ~any(tail)
            tail = m+1;
        end
        DG(head, tail) = 1;
        if rev(i)
            DG(tail, head) = 1;
        end
    end
    DG = sparse(DG);
    [~, C] = graphconncomp(DG, 'Directed', true);
    C = C(1:end-1);
    if range(C) == 0
        [~, C] = graphconncomp(max(DG, DG.'), 'Directed', false);
        C = C(1:end-1);
    end
    component = zeros(n, 1);
    for i = 1:n
        v = C(S(:, i) ~= 0);
        component(i) = all(v == v(1))*v(1);
    end
    if range(component) == 0
        component = zeros(n, 1);
        if strcmp(algorithm, 'fast')
            component(fastcc(model, 1e-4)) = 1;
        elseif strcmp(algorithm, 'swift')
            component(swiftcc(S, rev, solver)) = 1;
        end
    else
        newcomponent = zeros(n, 1);
        for i = unique(component).'
            core = component == i;
            if sum(core) > 1
                model.S = S(:, core);
                model.rev = rev(core);
                model.lb = lb(core);
                model.ub = ub(core);
                model.c = c(core);
                model.rxns = rxns(core);
                newcomponent(core) = partition(model, solver, algorithm);
            end
        end
        component = newcomponent;
    end
end