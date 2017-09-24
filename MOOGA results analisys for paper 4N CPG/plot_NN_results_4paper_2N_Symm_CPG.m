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

%% train NN 5 times and collect statistics about the Perf:
close all; clc;

% % get the names of the training parameters: (inverse)
% Inputs_names = {'period','tau','a'};
% Targets_names = {'b'};

% get the names of the training parameters: (farward)
Inputs_names = {'tau','b','a'};
Targets_names = {'period'};

architecture = {[20]};
numOfRepeats = 5;

train_RMSE = zeros(numOfRepeats,length(architecture));
valid_RMSE = zeros(numOfRepeats,length(architecture));
test_RMSE = zeros(numOfRepeats,length(architecture));

train_R2 = zeros(numOfRepeats,length(architecture));
valid_R2 = zeros(numOfRepeats,length(architecture));
test_R2 = zeros(numOfRepeats,length(architecture));

train_slope = zeros(numOfRepeats,length(architecture));
valid_slope = zeros(numOfRepeats,length(architecture));
test_slope = zeros(numOfRepeats,length(architecture));

for i=1:length(architecture)
    disp(['at NN arch :  ',num2str(architecture{1,i})]);
    for j=1:5
        disp(['    trail num #',num2str(j)])
        NNs_4paper = NNs_4paper.train_and_test(Inputs_names,Targets_names,...
        architecture{1,i},'NN',0);

        train_RMSE(j,i) = NNs_4paper.NN.train_RMSE;
        valid_RMSE(j,i) = NNs_4paper.NN.valid_RMSE;
        test_RMSE(j,i) = NNs_4paper.NN.test_RMSE;

        train_R2(j,i) = NNs_4paper.NN.train_R2;
        valid_R2(j,i) = NNs_4paper.NN.valid_R2;
        test_R2(j,i) = NNs_4paper.NN.test_R2;

        train_slope(j,i) = NNs_4paper.NN.train_slope;
        valid_slope(j,i) = NNs_4paper.NN.valid_slope;
        test_slope(j,i) = NNs_4paper.NN.test_slope;
    end
end

train_RMSE_mean = mean(train_RMSE,1);
valid_RMSE_mean = mean(valid_RMSE,1);
test_RMSE_mean = mean(test_RMSE,1)

train_R2_mean = mean(train_R2,1);
valid_R2_mean = mean(valid_R2,1);
test_R2_mean = mean(test_R2,1)

train_slope_mean = mean(train_slope,1);
valid_slope_mean = mean(valid_slope,1);
test_slope_mean = mean(test_slope,1)

figure;
boxplot(train_RMSE,'Colors',[0,0,128]./256); hold on;
boxplot(valid_RMSE,'Colors',[34,139,34]./256);
boxplot(test_RMSE,'Colors',[178,34,34]./256);
xlabel('NN arcith');
ylabel('RMSE');
grid minor
title('RMSE over hidden neurons num');

figure;
boxplot(train_R2,'Colors',[0,0,128]./256); hold on;
boxplot(valid_R2,'Colors',[34,139,34]./256);
boxplot(test_R2,'Colors',[178,34,34]./256);
xlabel('NN arcith');
ylabel('R^2');
grid minor
title('R^2 over hidden neurons num');

figure;
boxplot(train_slope,'Colors',[0,0,128]./256); hold on;
boxplot(valid_slope,'Colors',[34,139,34]./256);
boxplot(test_slope,'Colors',[178,34,34]./256);
xlabel('NN arcith');
ylabel('reggresion graph slope');
grid minor
title('slope over hidden neurons num');

