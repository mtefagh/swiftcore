clear all; close all; clc;
% importing the metabolic network model
load('consistRecon1.mat');
n = length(model.rxns);
solver = 'gurobi';
changeCobraSolver(solver);
% extracting the set of core reactions
core = load('C_liver.mat').C;
% setting the number of iterations
iterations = 10;
LP = zeros(2, iterations);
runtime = zeros(3, iterations);
performance = zeros(4, iterations);
errors = true(3, iterations);
finalScore = zeros(3, iterations);
for k = 1:iterations
    disp(k);
    % FASTCORE
    time1 = tic;
    coreRxn = fastcore(core, model, 1e-4);
    runtime(1, k) = toc(time1);
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
        tempmodel.rxnNames = model.rxnNames(coreRxnBool);
        tempmodel.rules = model.rules(coreRxnBool);
        tempmodel.c = model.c(coreRxnBool);
        tempmodel.mets = model.mets;
        tempmodel.b = model.b;
        A = fastcc(tempmodel, 1e-4);
        if all(A.' == 1:length(A))
            errors(1, k) = false;
        end
        taskReport = checkMetabolicTasks(tempmodel, 'liverTests.xlsx');
        taskReport = split(taskReport(211, 3), '/');
        finalScore(1, k) = str2double(taskReport(1));
    end
    % SWIFTCORE w/o reduction
    time2 = tic;
    [~, coreInd, numLP] = swiftcore(model, core, ones(n, 1), 1e-10, false, solver);
    runtime(2, k) = toc(time2);
    LP(1, k) = numLP;
    performance(2, k) = sum(coreInd);
    % consistency checking
    if all(coreInd(core))
        tempmodel.S = model.S(:, coreInd);
        tempmodel.rev = model.rev(coreInd);
        tempmodel.lb = model.lb(coreInd);
        tempmodel.ub = model.ub(coreInd);
        tempmodel.rxns = model.rxns(coreInd);
        tempmodel.rxnNames = model.rxnNames(coreInd);
        tempmodel.rules = model.rules(coreInd);
        tempmodel.c = model.c(coreInd);
        tempmodel.mets = model.mets;
        tempmodel.b = model.b;
        A = fastcc(tempmodel, 1e-4);
        if all(A.' == 1:length(A))
            errors(2, k) = false;
        end
        taskReport = checkMetabolicTasks(tempmodel, 'liverTests.xlsx');
        taskReport = split(taskReport(211, 3), '/');
        finalScore(2, k) = str2double(taskReport(1));
    end
    % SWIFTCORE w/ reduction
    time3 = tic;
    [~, coreIndReduce, numLP] = swiftcore(model, core, ones(n, 1), 1e-10, true, solver);
    runtime(3, k) = toc(time3);
    LP(2, k) = numLP;
    performance(3, k) = sum(coreIndReduce);
    % consistency checking
    if all(coreIndReduce(core))
        tempmodel.S = model.S(:, coreIndReduce);
        tempmodel.rev = model.rev(coreIndReduce);
        tempmodel.lb = model.lb(coreIndReduce);
        tempmodel.ub = model.ub(coreIndReduce);
        tempmodel.rxns = model.rxns(coreIndReduce);
        tempmodel.rxnNames = model.rxnNames(coreIndReduce);
        tempmodel.rules = model.rules(coreIndReduce);
        tempmodel.c = model.c(coreIndReduce);
        tempmodel.mets = model.mets;
        tempmodel.b = model.b;
        A = fastcc(tempmodel, 1e-4);
        if all(A.' == 1:length(A))
            errors(3, k) = false;
        end
        taskReport = checkMetabolicTasks(tempmodel, 'liverTests.xlsx');
        taskReport = split(taskReport(211, 3), '/');
        finalScore(3, k) = str2double(taskReport(1));
    end
    % calculating the size of the intersection of the answers of the three algorithms above
    performance(4, k) = sum(coreIndReduce + coreInd + coreRxnBool == 3);
end
% consistency checking
fprintf('FASTCORE returns wrong answers in %d%% of cases!\n', round(100*sum(errors(1,:))/size(errors,2)));
fprintf('SWIFTCORE w/o reduction returns wrong answers in %d%% of cases!\n', round(100*sum(errors(2,:))/size(errors,2)));
fprintf('SWIFTCORE w/ reduction returns wrong answers in %d%% of cases!\n', round(100*sum(errors(3,:))/size(errors,2)));
% drawing scatterplots for the comparison of the number of LPs
figure();
plot(1:iterations,LP(1,:),'o',1:iterations,LP(2,:),'*');
legend({'SWIFTCORE w/o reduction','SWIFTCORE w/ reduction'});
xlabel('iteration');
ylabel('number of solved LPs');
savefig('LPLiver');
% drawing scatterplots for the comparison of the runtime
figure();
plot(1:iterations,runtime(1,:),'o',1:iterations,runtime(2,:),'*',1:iterations,runtime(3,:),'+');
legend({'FASTCORE','SWIFTCORE w/o reduction','SWIFTCORE w/ reduction'},'Location','northeast');
xlabel('iteration');
ylabel('runtime in seconds');
savefig('runtimeLiver');
% drawing scatterplots for the comparison of the number of passed tests
figure();
plot(1:iterations,finalScore(1,:),'o',1:iterations,finalScore(2,:),'*',1:iterations,finalScore(3,:),'+');
legend({'FASTCORE','SWIFTCORE w/o reduction','SWIFTCORE w/ reduction'},'Location','northeast');
xlabel('iteration');
ylabel('number of passed tests');
savefig('testsLiver');
% drawing scatterplots for the comparison of the sparsity
figure();
plot(1:iterations,performance(1,:),'o',1:iterations,performance(2,:),'*',1:iterations,performance(3,:),'+',1:iterations,performance(4,:),'x');
legend({'FASTCORE','SWIFTCORE w/o reduction','SWIFTCORE w/ reduction','intersection'},'Location','southeast');
xlabel('iteration');
ylabel('size of the subnetwork');
savefig('performanceLiver');
% comparing the improvement percentage of the runtime
fprintf('switftcore w/o reduction was computed in %d%% of fastcore runtime.\n', round(mean(runtime(2, :)./runtime(1, :))*100));
fprintf('switftcore w/ reduction was computed in %d%% of fastcore runtime.\n', round(mean(runtime(3, :)./runtime(1, :))*100));
% comparing the improvement percentage of the sparsity
fprintf('switftcore w/o reduction performed with %d%% of fastcore reactions.\n', round(mean(performance(2, :)./performance(1, :))*100));
fprintf('switftcore w/ reduction performed with %d%% of fastcore reactions.\n', round(mean(performance(3, :)./performance(1, :))*100));
fprintf('The average number of reactions in the intersection is %d%% of fastcore reactions.\n', ...
    round(mean(performance(4, :)./performance(1, :))*100));