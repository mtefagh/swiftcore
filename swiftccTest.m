clear all; close all; clc;
% importing the metabolic network model
load('Recon3D_301/Recon3D_301.mat');
model = Recon3D;
clear Recon3D;
S = model.S;
rev = double(model.lb == -1000);
lb = model.lb*10000;
ub = model.ub*10000;
c = model.c;
rxns = model.rxns;
solver = 'cplex';
runtimecc = zeros(4, 29);
errors = false(1, 29);
for k = 1:29
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
    tic;
    A = fastcc(model, 1e-4);
    runtimecc(1, k) = toc;
    % SWIFTCC
    tic;
    consistent = swiftcc(model.S, model.rev, solver);
    runtimecc(2, k) = toc;
    % FASTCC++
    tic;
    component = partition(model, solver, 'fast');
    runtimecc(3, k) = toc;
    % SWIFTCC++
    tic;
    component2 = partition(model, solver, 'swift');
    runtimecc(4, k) = toc;
    % checking to see if all the four algorithms above return the same answer
    errors(1, k) = any(~component2(consistent)) || length(consistent) ~= sum(component2) ...
            || length(consistent) ~= length(A) || any(A ~= consistent) || any(component ~= component2);
end
if any(errors)
    warning('Wrong answer!');
end
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