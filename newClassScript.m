
clear all; close all; clc

fileName = 'MatsRandomRes_2Neurons_change_only_a.mat';
parametersCells = {'tau','b','a','s'};
targetCells = {'freq'};
seqOrder = {'tau','b','s','_','a'}; %'_' - dont care parameter
myCode = myCode(fileName,parametersCells,targetCells,seqOrder);

%% with NN:
myCode = myCode.Set('NN',[5,5]);
myCode = myCode.trainNN(1);

%% with my_MoE
myCode = myCode.Set('our_MoE',10,2,[2],[2],5,3);
myCode = myCode.my_MoE_train();
myCode.my_MoE_plot_train_perf();
close all;
myCode.my_MoE_testNet(myCode.sampl_test,myCode.targ_test,myCode.my_MoE_out.expertsNN,...
    myCode.my_MoE_out.gateNet,3,1);

% TODO: check tjis code for competetiveflag =1,2