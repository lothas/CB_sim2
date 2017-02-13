
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

%% 4N CPG symmetric
fileName = 'MatsRandomRes_4Neurons_symm.mat';
parametersCells = {'tau','b','w_12','w_13','w_14','w_23','w_24','w_34'};
targetCells = {'freq'};
seqOrder = {'tau','b','c','w_12','w_13','w_14','w_23','w_24','w_34'};
myCode = myCode(fileName,parametersCells,targetCells,seqOrder);

clear fileName parametersCells targetCells seqOrder

problemType = '4N_symm';
%% with NN:
myCode = myCode.Set('NN',[5],100);
myCode = myCode.trainNN(0);
myCode.plot_fit_data('NN',problemType);
%% with my_MoE
competetiveflag = 3;
myCode = myCode.Set('our_MoE',5,3,[2],[2],5,competetiveflag);
myCode = myCode.my_MoE_train();
myCode.my_MoE_plot_train_perf();
myCode.my_MoE_testNet(myCode.sampl_test,myCode.targ_test,myCode.my_MoE_out.expertsNN,...
    myCode.my_MoE_out.gateNet,competetiveflag,1);
myCode.plot_fit_data('my_MoE',problemType);

%% with the paper's MoE
myCode = myCode.Set('paper_MoE',40,5,0.005,0.995);
myCode = myCode.paper_MoE_train();
myCode.paper_MoE_plotPerf({'MSE_over_iter','regGraph','reg_with_color','reg_graph_from_NNtoolbox'})
myCode.plot_fit_data('paper_MoE',problemType);


