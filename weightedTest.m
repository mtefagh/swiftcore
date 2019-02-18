clear all; close all; clc;
% importing the metabolic network model
load('Recon3D_301/Recon3DModel_301.mat');
model = Recon3DModel;
model.rev = double(model.lb < 0);
clear Recon3DModel;
n = length(model.rxns);
% setting up the LP solver to IBM CPLEX
solver = 'ibm_cplex';
changeCobraSolver(solver);
% selecting a range of different magnitudes for the weights
weights = 2.^(-4:4);
l = length(weights);
LP = zeros(2, l);
runtime = zeros(3, l);
performance = zeros(4, l);
errors = true(3, l);
for k = 1:l
    disp(k);
    % randomly sampling the set of core reactions
    core = randsample(n, n/2);
    % FASTCORE
    tic;
    coreRxn = fastCoreWeighted(core, model, weights(k)*ones(n, 1), 1e-4);
    runtime(1, k) = toc;
    coreRxnBool = false(n, 1);
    coreRxnBool(coreRxn) = true;
    performance(1, k) = sum(coreRxnBool);
    % consistency checking
    if all(coreRxnBool(core))
        tempmodel.S = model.S(:, coreRxnBool);
        tempmodel.rev = model.rev(coreRxnBool);
        tempmodel.lb = model.lb(coreRxnBool);
        tempmodel.ub = model.ub(coreRxnBool);
        tempmodel.rxns = model.rxns(coreRxnBool);
        tempmodel.c = model.c(coreRxnBool);
        tempmodel.mets = model.mets;
        A = swiftcc(tempmodel.S, tempmodel.rev, solver);
        if all(A.' == 1:length(A))
            errors(1, k) = false;
        end
    end
    % SWIFTCORE w/o reduction
    tic;
    [coreInd, numLP] = swiftcore(model.S, model.rev, core, weights(k)*ones(n, 1), false, solver);
    runtime(2, k) = toc;
    LP(1, k) = numLP;
    performance(2, k) = sum(coreInd);
    % consistency checking
    if all(coreInd(core))
        tempmodel.S = model.S(:, coreInd);
        tempmodel.rev = model.rev(coreInd);
        tempmodel.lb = model.lb(coreInd);
        tempmodel.ub = model.ub(coreInd);
        tempmodel.rxns = model.rxns(coreInd);
        tempmodel.c = model.c(coreInd);
        tempmodel.mets = model.mets;
        A = swiftcc(tempmodel.S, tempmodel.rev, solver);
        if all(A.' == 1:length(A))
            errors(2, k) = false;
        end
    end
    % SWIFTCORE w/ reduction
    tic;
    [coreIndReduce, numLP] = swiftcore(model.S, model.rev, core, weights(k)*ones(n, 1), true, solver);
    runtime(3, k) = toc;
    LP(2, k) = numLP;
    performance(3, k) = sum(coreIndReduce);
    % consistency checking
    if all(coreIndReduce(core))
        tempmodel.S = model.S(:, coreIndReduce);
        tempmodel.rev = model.rev(coreIndReduce);
        tempmodel.lb = model.lb(coreIndReduce);
        tempmodel.ub = model.ub(coreIndReduce);
        tempmodel.rxns = model.rxns(coreIndReduce);
        tempmodel.c = model.c(coreIndReduce);
        tempmodel.mets = model.mets;
        A = swiftcc(tempmodel.S, tempmodel.rev, solver);
        if all(A.' == 1:length(A))
            errors(3, k) = false;
        end
    end
    % calculating the size of the intersection of the answers of the three algorithms above
    performance(4, k) = sum(coreIndReduce + coreInd + coreRxnBool == 3);
end
% consistency checking
if any(errors, 'all')
    warning('Wrong answers!');
end
% drawing boxplots for the comparison of runtimes
figure();
boxplot(runtime.','Notch','on','Labels',{'FASTCORE','SWIFTCORE w/o reduction','SWIFTCORE w/ reduction'});
ylabel('runtime in seconds');
savefig('runtimeBox');
% drawing boxplots for the comparison of sparsity
figure();
boxplot(performance.','Notch','on','Labels',{'FASTCORE','SWIFTCORE w/o reduction','SWIFTCORE w/ reduction','intersection'});
ylabel('size of the subnetwork');
savefig('performanceBox');
% drawing boxplots for the comparison of the number of LPs
figure();
boxplot(LP.','Notch','on','Labels',{'SWIFTCORE w/o reduction','SWIFTCORE w/ reduction'});
ylabel('number of solved LPs');
savefig('LPBox');