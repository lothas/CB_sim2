
%% NN results:
clear all; close all; clc

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 2;
% % change tau_a/tau_r to 12 (instead of 5)
MML.Sim.Con.tau_ratio = 12;

seqOrder = {'tau','b','c','NR','a','NR','NR','NR'}; % 'NR'-not relevant

% Load data:
% results_fileName = 'MatsRandomRes_2Neurons_symm_4Paper_All_in_range.mat';
results_fileName = 'MatsRandomRes_4Neurons_Large_b_Large_W_All_osc';
load(results_fileName,'results');

% get CPG periods:
periods = horzcat(results(:).periods);

% get the sequences:
seq = (vertcat(results(:).seq))';
%% Stage#1: filter results
osc_ids_temp = ~isnan(periods);
disp(['Number of non-osc CPGs: ',num2str(sum(~osc_ids_temp))]);

ids_in_genome_range = true(1,size(seq,2));
for n=1:MML.Gen.Length
    if ~strcmp(seqOrder{1,n},'NR')
        ids_temp = (seq(n,:) > MML.Gen.Range(1,n)) &...
            (seq(n,:) < MML.Gen.Range(2,n));
        ids_in_genome_range = ids_in_genome_range & ids_temp;
    end
end

disp(['Number of CPGs with parameters not in range: ',...
    num2str(sum(~ids_in_genome_range))]);

osc_ids = osc_ids_temp & ids_in_genome_range;

% get CPG ids that are in Period range:
osc_inRange_ids = osc_ids &...
    ( (periods(1,:) > MML.perLimOut(1,1)) &...
    (periods(1,:) < MML.perLimOut(1,2)) );

%% stage 2: make NN in and out:
inputs = [periods;...
    seq(1,:);...
    seq(5,:)];

targets = seq(2,:);


%% Stage 3a: use "normal" NN
% Create and train the NN
net = feedforwardnet(5);

% % specilize train functions:
% net.trainFcn = 'trainbr';
% net.trainFcn = 'trainscg';
% net.trainFcn = 'trainrp';

[net, ~] = train(net, inputs, targets);

outputs = net(inputs);

figure;
histogram2(targets,outputs,100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
xlabel('targets');
ylabel('outputs');
axis equal

figure;
scatter(targets,outputs);
    
%% Stage 3b: use Matlab's other toolbox:

% define NN layers
layers = [reluLayer
    fullyConnectedLayer(1)
    regressionLayer];

% training options:
options = trainingOptions('sgdm','InitialLearnRate',0.001, ...
    'MaxEpochs',15);


net = trainNetwork(inputs,targets,layers,options);

outputs = predict(net,inputs);

% calc RMSE:
err = targets - outputs;
sqrErr = err .^ 2;
RMSE = sqrt(mean(sqrErr))

figure;
histogram2(targets,outputs,100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
xlabel('targets');
ylabel('outputs');
axis equal

figure;
scatter(targets,outputs);