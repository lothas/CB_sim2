%% NN results:
clear all; close all; clc

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;

% % % data with many CPG's that oscillates in range:
% results_fileName = 'MatsRandomRes_4Neurons_4Paper_for_MOOGA_try.mat';

% % data with small amount of CPGs that oscillate in range:
% results_fileName = 'MatsRandomRes_4Neurons_4Paper.mat';

% % data with tau_ratio=12 and (0.2 < b < 2.5)
results_fileName = 'MatsRandomRes_4Neurons_4Paper_tau_ratio_equalTo_12_added_b_4Paper1.mat';
% % change tau_a/tau_r to 12 (instead of 5)
MML.Sim.Con.tau_ratio = 12;
% % change b_max:
MML.Gen.Range(2,2) = 2.5; % the class will filter genes that are not in the new range.

NNs_4paper = NNs_4paper(results_fileName,MML);

%% plot histograms of  osc_param Vs n-osc_param:
paramName = 'tau';
norm_flag = true;
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
architecture = [20,20];
NNs_4paper.train_and_test_NN(caseNum,architecture,'test_on_training_data');

%% correlation between targets and NN outputs:
% a = [tau_out_on_train_set;b_out_on_train_set];
% b = [tau_training;b_training];
% [RHO,PVAL] = corr(a',b')