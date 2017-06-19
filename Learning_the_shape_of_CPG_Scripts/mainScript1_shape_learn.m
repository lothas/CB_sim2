
clc; clear all; close all;

%% define the class:
filename = 'MatsRandomRes_2Neurons_symm_test_for_FFT_only_conv.mat';
% inputsNames = {'tau','b','a','c'};
% outputNames = {'A0','A1','A2','A3','B1','B2','B3','freq'};

inputsNames = {'c','tau'};
outputNames = {'A1','B1'};

CPGshapeLearn = CPGshapeLearn(filename,inputsNames,outputNames);

%% NN training:
hidN = [20,20];
epochN = 500;

CPGshapeLearn = CPGshapeLearn.NN_training(hidN,epochN,1,0);

errs = CPGshapeLearn.NN.MSE_test_perf;
disp('the MSE err for each output:');
for i=1:length(outputNames)
    disp([outputNames{1,i},' = ',num2str(errs(i))]);
end

CPGshapeLearn.plot_param_dis('all')
% TODO: add training for each targ seperatly
