clear all; close all; clc;
% importing the metabolic network model
load('Recon3D_301/Recon3DModel_301.mat');
model = Recon3DModel;
model.rev = double(model.lb < 0);
clear Recon3DModel;
n = length(model.rxns);
solver = 'ibm_cplex';
LP = zeros(2, 105);
runtime = zeros(3, 105);
performance = zeros(4, 105);
errors = true(3, 105);
for k = 1:105
    disp(k);
    % randomly sampling the set of core reactions
    core = randsample(n, 100*k);
    % FASTCORE
    tic;
    coreRxn = fastcore(core, model, 1e-4);
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
    [coreInd, numLP] = swiftcore(model.S, model.rev, core, ones(n, 1), false, solver);
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
    [coreIndReduce, numLP] = swiftcore(model.S, model.rev, core, ones(n, 1), true, solver);
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
fprintf('FASTCORE returns wrong answers in %d%% of cases!\n', round(100*sum(errors(1,:))/size(errors,2)));
fprintf('SWIFTCORE w/o reduction returns wrong answers in %d%% of cases!\n', round(100*sum(errors(2,:))/size(errors,2)));
fprintf('SWIFTCORE w/ reduction returns wrong answers in %d%% of cases!\n', round(100*sum(errors(3,:))/size(errors,2)));
% drawing scatterplots for the comparison of the number of LPs
figure();
plot(100*(1:105),LP(1,:),'o',100*(1:105),LP(2,:),'*');
legend({'SWIFTCORE w/o reduction','SWIFTCORE w/ reduction'});
xlabel('number of core reactions');
ylabel('number of solved LPs');
savefig('LP');
% drawing scatterplots for the comparison of the runtime
figure();
plot(100*(1:105),runtime(1,:),'o',100*(1:105),runtime(2,:),'*',100*(1:105),runtime(3,:),'+');
legend({'FASTCORE','SWIFTCORE w/o reduction','SWIFTCORE w/ reduction'},'Location','northeast');
xlabel('number of core reactions');
ylabel('runtime in seconds');
savefig('runtime');
% drawing scatterplots for the comparison of the sparsity
figure();
plot(100*(1:105),performance(1,:),'o',100*(1:105),performance(2,:),'*',...
    100*(1:105),performance(3,:),'+',100*(1:105),performance(4,:),'x');
legend({'FASTCORE','SWIFTCORE w/o reduction','SWIFTCORE w/ reduction','intersection'},'Location','southeast');
xlabel('number of core reactions');
ylabel('size of the subnetwork');
savefig('performance');
% comparing the improvement percentage of the runtime
fprintf('switftcore w/o reduction was computed in %d%% of fastcore runtime.\n', round(mean(runtime(2, :)./runtime(1, :))*100));
fprintf('switftcore w/ reduction was computed in %d%% of fastcore runtime.\n', round(mean(runtime(3, :)./runtime(1, :))*100));
% comparing the improvement percentage of the sparsity
fprintf('switftcore w/o reduction performed with %d%% of fastcore reactions.\n', round(mean(performance(2, :)./performance(1, :))*100));
fprintf('switftcore w/ reduction performed with %d%% of fastcore reactions.\n', round(mean(performance(3, :)./performance(1, :))*100));
fprintf('The average number of reactions in the intersection is %d%% of fastcore reactions.\n', ...
    round(mean(performance(4, :)./performance(1, :))*100));