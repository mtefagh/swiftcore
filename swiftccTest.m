clear all; close all; clc;
% importing the metabolic network model
changeCobraSolver('ibm_cplex');
load('Recon3D_301/Recon3D_301.mat');
model = Recon3D;
clear Recon3D;
S = model.S;
rev = double(model.lb == -1000);
lb = model.lb;
ub = model.ub;
c = model.c;
rxns = model.rxns;
solver = 'cplex';
runtimecc = zeros(4, 29);
errors = false(2, 29);
unknown_error = false;
for k = 1:29
    error = false;
    disp(k);
    % randomly selecting a sample subnetwork of the given size
    core = randsample(length(rxns), 467*k);
    model.S = S(:, core);
    model.rev = rev(core);
    model.lb = lb(core);
    model.ub = ub(core);
    model.c = c(core);
    model.rxns = rxns(core);
    % FASTCC
    time1 = tic;
    A = fastcc(model, 1e-4);
    runtimecc(1, k) = toc(time1);
    % SWIFTCC
    time2 = tic;
    consistent = swiftcc(model.S, model.rev, solver);
    runtimecc(2, k) = toc(time2);
    % consistency checking
    if length(A) ~= length(consistent)|| any(A ~= consistent)
        error = true;
    end
    % FASTCC++
    time3 = tic;
    component = partition(model, solver, 'fast');
    runtimecc(3, k) = toc(time3);
    % consistency checking
    if any(~component(consistent)) || length(consistent) ~= sum(component)
        error = true;
    end
    % SWIFTCC++
    time4 = tic;
    component2 = partition(model, solver, 'swift');
    runtimecc(4, k) = toc(time4);
    % consistency checking
    if length(component) ~= length(component2) || any(component ~= component2)
        error = true;
    end
    % checking to see if all the four algorithms above return the same answer
    if error
        if ~isempty(identifyFastBlockedRxns(model, model.rxns(A)))
            errors(1, k) = true;
        elseif ~isempty(identifyFastBlockedRxns(model, model.rxns(consistent)))
            errors(2, k) = true;
        else
            unknown_error = true;
        end
    end
end
% consistency checking
if unknown_error
    warning('Wrong answers!');
end
fprintf('FASTCC returns wrong answers in %d%% of cases!\n', round(100*sum(errors(1,:))/size(errors,2)));
fprintf('SWIFTCC returns wrong answers in %d%% of cases!\n', round(100*sum(errors(2,:))/size(errors,2)));
% drawing scatterplots for the comparison of runtimes
figure();
plot(467*(1:29),runtimecc(1,:),'o',467*(1:29),runtimecc(2,:),'*',467*(1:29), ...
    runtimecc(3,:),'+',467*(1:29),runtimecc(4,:),'x');
legend({'FASTCC','SWIFTCC','FASTCC++','SWIFTCC++'},'Location','northwest');
xlabel('number of core reactions');
ylabel('runtime in seconds');
savefig('runtimeccScatter');
% drawing boxplots for the comparison of runtimes
figure();
boxplot(runtimecc.','Notch','on','Labels',{'FASTCC','SWIFTCC','FASTCC++','SWIFTCC++'});
ylabel('Runtime');
savefig('runtimeccBox');
% comparing the improvement percentage of the runtime
fprintf('switftcc was computed in %d%% of fastcc runtime.\n', round(sum(runtimecc(2, :))/sum(runtimecc(1, :))*100));
fprintf('fastcc++ was computed in %d%% of fastcc runtime.\n', round(sum(runtimecc(3, :))/sum(runtimecc(1, :))*100));
fprintf('switftcc++ was computed in %d%% of fastcc runtime.\n', round(sum(runtimecc(4, :))/sum(runtimecc(1, :))*100));