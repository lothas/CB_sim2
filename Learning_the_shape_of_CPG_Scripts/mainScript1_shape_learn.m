
clc; clear all; close all;

%% define the class:
filename = 'MatsRandomRes_2Neurons_symm_test_for_FFT_only_conv.mat';
% inputsNames = {'tau','b','a','c'};
% outputNames = {'A0','A1','A2','A3','B1','B2','B3','freq'};

inputsNames = {'c','b','a'};
outputNames = {'A1','B1'};

CPGshapeLearn = CPGshapeLearn(filename,inputsNames,outputNames);

%% NN training:
clc 

hidN = 10;
epochN = 500;

CPGshapeLearn = CPGshapeLearn.NN_training(hidN,epochN,1,0);

%% multi NN training:
clc 

hidN = 50;
epochN = 500;
CPGshapeLearn = CPGshapeLearn.NN_multi_training(hidN,epochN,0);

%%
CPGshapeLearn.plot_param_dis('all')

CPGshapeLearn.plot_param_dis('B1')
% TODO: add training for each targ seperatly
