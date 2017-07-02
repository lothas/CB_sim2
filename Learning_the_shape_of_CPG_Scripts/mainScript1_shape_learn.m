
clc; clear all; close all;

%% define the class:
% filename = 'MatsRandomRes_2Neurons_symm_test_for_FFT_only_conv.mat';
filename = 'MatsRandomRes_2Neurons_symm_test_for_shapeLearn_changing_C.mat';
% 
% % % % Farward
% inputsNames = {'tau','b','a','c'};
% outputNames = {'A0','A1','A2','A3','B1','B2','B3','freq'};

% % % % Reverse
% inputsNames = {'A0','A1','A2','A3','B1','B2','B3','freq'};
% outputNames = {'tau','b','a','c'};

% % Custom:
inputsNames = {'A1','B1','b','a','tau','freq'};
outputNames = {'c'};

CPGshapeLearn = CPGshapeLearn(filename,inputsNames,outputNames);

%% NN training:
clc 

hidN = 10;
epochN = 500;

CPGshapeLearn = CPGshapeLearn.NN_training(hidN,epochN,1,0);

test_ind = CPGshapeLearn.NN.net_perf.testMask{1,1};
test_ind = ~isnan(test_ind);
NN_out = CPGshapeLearn.NN.net(CPGshapeLearn.sampl);
NN_test_out = NN_out(test_ind);
NN_test_targ = CPGshapeLearn.targ(test_ind);
plotregression(NN_test_targ, NN_test_out, 'Testing')
%% multi NN training:
clc 

hidN = [10,10,10,10];
epochN = 1000;
CPGshapeLearn = CPGshapeLearn.NN_multi_training(hidN,epochN,0);

% NNs = CPGshapeLearn.multi_NN.NNs;
% NNs_tr = CPGshapeLearn.multi_NN.net_perf;
% outNames = CPGshapeLearn.multi_NN.outNames;
% 
% num_of_NNs = length(NNs);
% 
% j = 8;
% 
% % plotperform(NNs_tr{1,j})
% 
% NN = NNs{1,j};
% tr = NNs_tr{1,j};
% test_ind = tr.testMask{1,1};
% test_ind = ~isnan(test_ind);
% NN_out = NN(CPGshapeLearn.sampl);
% NN_test_out = NN_out(test_ind);
% NN_test_targ = CPGshapeLearn.targ(test_ind);
% plotregression(NN_test_targ, NN_test_out, 'Testing')
%%
CPGshapeLearn.plot_param_dis('all')

CPGshapeLearn.plot_param_dis('A1')
% TODO: add training for each targ seperatly
