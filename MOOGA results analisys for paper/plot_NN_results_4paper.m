%% NN results:
clear all; close all; clc

% data with many CPG's that oscillates in range:
% results_fileName = 'MatsRandomRes_4Neurons_4Paper_for_MOOGA_try.mat';

% data with small amount of CPGs that oscillate in range:
results_fileName = 'MatsRandomRes_4Neurons_4Paper.mat';

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;

NNs_4paper = NNs_4paper(results_fileName,MML);

%% plot histograms of  osc_param Vs n-osc_param:
paramName = 'tau';
NNs_4paper.plot_oscParam_vs_NoscParam_hist(paramName)

%% rescale tau
NNs_4paper = NNs_4paper.rescale_CPGs();

tau_before = NNs_4paper.seq(:,1);
tau_after = NNs_4paper.tau_rescaled;

% norm results with min,max:
tau_before = NNs_4paper.norm_min_max(tau_before,'tau');
tau_after = NNs_4paper.norm_min_max(tau_after,'tau');

NNs_4paper.hist_compare(tau_before,tau_after,'tau',...
    20,{'tau','tau_{rescaled}'},'plot')

%% train NN:

% make train groups names:
caseNum = 9;
[Inputs_names,Targets_names] = NNs_4paper.check_NN_case(caseNum,'period');
[Inputs,Targets] = NNs_4paper.prepare_NN_train_data(Inputs_names,Targets_names);

architecture = [10];
NNs_4paper = NNs_4paper.train_NN(architecture,Inputs,Targets);

[Inputs_names,Targets_names] = ...
    NNs_4paper.check_NN_case(caseNum,'period_desired');

seq_test  = MML.Gen.RandSeq(1000);

seq_after_NN = ...
    NNs_4paper.prepare_NN_test_seq(seq_test,Inputs_names,Targets_names);

% compare tau before or after NN:
tau_before = seq_test(:,1);
tau_after = seq_after_NN(:,1);

% norm results with min,max:
tau_before = NNs_4paper.norm_min_max(tau_before,'tau');
tau_after = NNs_4paper.norm_min_max(tau_after,'tau');

NNs_4paper.hist_compare(tau_before,tau_after,'tau',...
    20,{'tau','tau_{rescaled}'},'plot');



