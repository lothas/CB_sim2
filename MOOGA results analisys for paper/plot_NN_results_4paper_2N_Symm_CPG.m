%% NN results:
clear all; close all; clc

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 2;

% % data with tau_ratio=12 and (0.2 < b < 2.5)
results_fileName = 'MatsRandomRes_2Neurons_symm_4Paper.mat';
% results_fileName = 'MatsRandomRes_2Neurons_symm_osc_in_range.mat';
% results_fileName = 'MatsRandomRes_2Neurons_symm_4Paper_All_in_range.mat';
% % change tau_a/tau_r to 12 (instead of 5)
MML.Sim.Con.tau_ratio = 12;

NNs_4paper = NNs_4paper(results_fileName,MML);

%% plot histograms of  osc_param Vs n-osc_param:
paramName = 'tau';
norm_flag = false;
NNs_4paper.plot_oscParam_vs_NoscParam_hist(paramName,norm_flag)
NNs_4paper.plot_oscParam_vs_oscInRangeParam_hist(paramName,norm_flag);

paramName = 'b';
norm_flag = false;
NNs_4paper.plot_oscParam_vs_NoscParam_hist(paramName,norm_flag);
NNs_4paper.plot_oscParam_vs_oscInRangeParam_hist(paramName,norm_flag);
%% plot 2D histogram:
norm_flag = true;
figure;
NNs_4paper.plot_2D_hist({'tau','b'},norm_flag)

periods = horzcat(NNs_4paper.results(NNs_4paper.osc_ids).periods);
figure;
histogram(periods,100,'Normalization','pdf');
xlabel('period [sec]');
title('histogram of the CPGs period');

%% rescale tau
NNs_4paper = NNs_4paper.rescale_CPGs();

tau_before = NNs_4paper.seq(:,1);
tau_after = NNs_4paper.tau_rescaled;

% norm results with min,max:
tau_before = NNs_4paper.norm_min_max(tau_before,'tau');
tau_after = NNs_4paper.norm_min_max(tau_after,'tau');

NNs_4paper.hist_compare(tau_before,tau_after,'tau',...
    20,{'tau','tau_{rescaled}'},'plot')

%% check NN on training data:
close all; clc;

caseNum = 7;

% get the names of the training parameters:
Inputs_names = {'period','tau','a'};
Targets_names = {'b'};

architecture = [20,20];

NNs_4paper = NNs_4paper.train_and_test(Inputs_names,Targets_names,...
    architecture,'NN','test_on_training_data');

% NNs_4paper.train_and_test(Inputs_names,Targets_names,...
%     architecture,'MoE colaboration','test_on_training_data');
% 
% NNs_4paper.train_and_test(Inputs_names,Targets_names,...
%     architecture,'MoE hard','test_on_training_data');
% 
% NNs_4paper.train_and_test(Inputs_names,Targets_names,...
%     architecture,'MoE soft','test_on_training_data');

%% save data to CSV file:
caseNum = 9;
fileName = 'data_to_CSV_case_7_tauRatio_12_';
NNs_4paper.save_data_CSV(caseNum,fileName);
%% correlation between targets and NN outputs:
% a = [tau_out_on_train_set;b_out_on_train_set];
% b = [tau_training;b_training];
% [RHO,PVAL] = corr(a',b')