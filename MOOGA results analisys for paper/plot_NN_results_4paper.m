%% NN results:
clear all; close all; clc

% data with many CPG's that oscillates in range:
% results_fileName = 'MatsRandomRes_4Neurons_4Paper_for_MOOGA_try.mat';

% data with small amount of CPGs that oscillate in range:
results_fileName = 'MatsRandomRes_4Neurons_4Paper.mat';
% results_fileName = 'MatsRandomRes_4Neurons_4Paper_for_MOOGA_try.mat';

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
norm_flag = true;

NNs_4paper.plot_oscParam_vs_NoscParam_hist(paramName,norm_flag)

NNs_4paper.plot_oscParam_vs_oscInRangeParam_hist(paramName,norm_flag)
%% plot 2D histogram:
norm_flag = true;
NNs_4paper.plot_2D_hist({'tau','b'},norm_flag)


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

% get param in taining data:
tau_training = Targets(1,:); % only for case #9 !
b_training = Targets(2,:); % only for case #9 !

architecture = [15];
NNs_4paper = NNs_4paper.train_NN(architecture,Inputs,Targets);

[Inputs_names,Targets_names] = ...
    NNs_4paper.check_NN_case(caseNum,'period_desired');

seq_test  = MML.Gen.RandSeq(100000);

seq_after_NN = ...
    NNs_4paper.prepare_NN_test_seq(seq_test,Inputs_names,Targets_names);

% compare tau before or after NN:
tau_before = seq_test(:,1);
tau_after = seq_after_NN(:,1);
b_before = seq_test(:,2);
b_after = seq_after_NN(:,2);

% norm results with min,max:
tau_training = NNs_4paper.norm_min_max(tau_training,'tau');
tau_before = NNs_4paper.norm_min_max(tau_before,'tau');
tau_after = NNs_4paper.norm_min_max(tau_after,'tau');
b_training = NNs_4paper.norm_min_max(b_training,'b');
b_before = NNs_4paper.norm_min_max(b_before,'b');
b_after = NNs_4paper.norm_min_max(b_after,'b');

binsNum = 20;


NNs_4paper.hist_compare(tau_before,tau_after,'tau',...
    binsNum,{'tau','tau_{after NN}'},'plot');

NNs_4paper.hist_compare(b_before,b_after,'b',...
    binsNum,{'b','b_{after NN}'},'plot');

% % plot 2D hist: (before NN)
% figure;
% histogram2(tau_before,b_before,binsNum,...
%     'DisplayStyle','tile','ShowEmptyBins','on',...
%     'Normalization','pdf');
% title('2D distribution of \tau and b before the NN');
% xlabel('tau norm');
% ylabel('b norm');
% axis([0,1,0,1]);
% 
% % plot 2D hist: (after NN)figure;
% histogram2(tau_after,b_after,binsNum,...
%     'DisplayStyle','tile','ShowEmptyBins','on',...
%     'Normalization','pdf');
% title('2D distribution of \tau and b after the NN');
% xlabel('tau norm');
% ylabel('b norm');
% axis([0,1,0,1]);

Xedges = linspace(0,1,100);
Yedges = linspace(0,1,100);
% plot 2D hist: (before NN)
figure;
histogram2(tau_before,b_before,Xedges,Yedges,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of tau and b before the NN');
xlabel('tau norm');
ylabel('b norm');
axis([0,1,0,1]);

% plot 2D hist: (after NN)figure;
histogram2(tau_after,b_after,Xedges,Yedges,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of tau and b after the NN');
xlabel('tau norm');
ylabel('b norm');
axis([0,1,0,1]);

% plot 2D hist: (after NN)figure;
figure; hold on
histogram2(tau_training,b_training,Xedges,Yedges,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
% scatter(tau_training,b_training,'.k');
title(sprintf('2D distribution of tau and b  \n in training data'));
xlabel('tau norm');
ylabel('b norm');
axis([0,1,0,1]);

figure;
histogram(tau_training,100);

figure;
histogram(b_training,100);

%% checking NN on the same training data set:

NNoutputs_on_train_set = NNs_4paper.NN.net(Inputs);
tau_out_on_train_set = NNs_4paper.norm_min_max(NNoutputs_on_train_set(1,:),'tau');
b_out_on_train_set = NNs_4paper.norm_min_max(NNoutputs_on_train_set(2,:),'b');

figure; hold on;
% scatter(tau_out_on_train_set,b_out_on_train_set,'.k');
histogram2(tau_out_on_train_set,b_out_on_train_set,Xedges,Yedges,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
scatter(tau_out_on_train_set,b_out_on_train_set,'.k');
title(sprintf('2D distribution of tau and b  \n on training data'));
xlabel('tau norm');
ylabel('b norm');
axis([0,1,0,1]);

%% correlation between targets and NN outputs:
a = [tau_out_on_train_set;b_out_on_train_set];
b = [tau_training;tau_training];
[RHO,PVAL] = corr(a',b')