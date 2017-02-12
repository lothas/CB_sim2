
clear all; close all; clc


%% 2N symmetric
fileName = 'MatsRandomRes_2Neurons_change_only_a.mat';
parametersCells = {'tau','b','a','s'};
targetCells = {'freq'};
seqOrder = {'tau','b','s','_','a'}; %'_' - dont care parameter
myCode = myCode(fileName,parametersCells,targetCells,seqOrder);

clear fileName parametersCells targetCells seqOrder

problemType = '2N_symm';
%% 2N general
fileName = 'MatsRandomRes_2Neurons_general_1.mat';
parametersCells = {'tau','b','w_12','w_21'};
targetCells = {'freq'};
seqOrder = {'tau','b','c_1','c_2','w_12','w_21'};
myCode = myCode(fileName,parametersCells,targetCells,seqOrder);

clear fileName parametersCells targetCells seqOrder

problemType = '2N_general';
%% with NN:
myCode = myCode.Set('NN',[5],100);
myCode = myCode.trainNN(1);

%% with my_MoE
competetiveflag = 2;
myCode = myCode.Set('our_MoE',10,2,[2],[2],5,competetiveflag);
myCode = myCode.my_MoE_train();
myCode.my_MoE_plot_train_perf();
close all;
myCode.my_MoE_testNet(myCode.sampl_test,myCode.targ_test,myCode.my_MoE_out.expertsNN,...
    myCode.my_MoE_out.gateNet,competetiveflag,1);

%% with the paper's MoE
myCode = myCode.Set('paper_MoE',10,2,0.001,0.995);
myCode = myCode.paper_MoE_train();
myCode.paper_MoE_plotPerf({'MSE_over_iter','regGraph','reg_graph_from_NNtoolbox'})

%%
out_net_test = myCode.NN.out_from_test;
% out_net_test = myCode.my_MoE_out.out_from_test;
% out_net_test = myCode.paper_MoE_out.out_from_test;
myCode.plot_fit_data(out_net_test,problemType);