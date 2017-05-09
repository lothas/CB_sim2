
clear all; close all; clc

%% 2N symmetric
% fileName = 'MatsRandomRes_2Neurons_change_only_a.mat';
% fileName = {'MatsRandomRes_2Neurons_symm_trainData_wide_range.mat',...
%     'MatsRandomRes_2Neurons_symm_testData.mat'};
fileName = {'MatsRandomRes_2Neurons_symm_trainData_wide_range_all.mat',...
    'MatsRandomRes_2Neurons_symm_testData2.mat'};
% fileName = {'MatsRandomRes_2Neurons_symm_trainData_narrow_range.mat',...
%     'MatsRandomRes_2Neurons_symm_testData.mat'};
parametersCells = {'tau','b','a','s'};
targetCells = {'freq'};
seqOrder = {'tau','b','s','_','a'}; %'_' - dont care parameter
myCode = myCode(fileName,parametersCells,targetCells,seqOrder,2);
myCode.disp_information = true;
myCode.sizeOfCPG = 2;
% clear fileName parametersCells targetCells seqOrder

problemType = '2N_symm';
%% 2N general
fileName = 'MatsRandomRes_2Neurons_general_1to4_combained.mat';
parametersCells = {'tau','b','w_{12}','w_{21}'};
targetCells = {'freq'};
seqOrder = {'tau','b','c_1','c_2','w_{12}','w_{21}'};
myCode = myCode(fileName,parametersCells,targetCells,seqOrder,2,true);
myCode.sizeOfCPG = 2;
% clear fileName parametersCells targetCells seqOrder

problemType = '2N_general';

%% 4N CPG symmetric
fileName = 'MatsRandomRes_4Neurons_symm_all_samples.mat';
% fileName = 'MatsRandomRes_4Neurons_symm.mat';
parametersCells = {'tau','b','w_{12}','w_{13}','w_{14}','w_{23}','w_{24}','w_{34}'};
targetCells = {'freq'};
seqOrder = {'tau','b','c','w_{12}','w_{13}','w_{14}','w_{23}','w_{24}','w_{34}'};
myCode = myCode(fileName,parametersCells,targetCells,seqOrder,4);
myCode.sizeOfCPG = 4;
% clear fileName parametersCells targetCells seqOrder

problemType = '4N_symm';

%% 4N CPG general
% fileName = 'MatsRandomRes_test.mat';
% fileName = 'MatsRandomRes_31_01_2017.mat';
fileName = 'MatsRandomRes_2-4_02_2017.mat';
parametersCells = {'tau','b',...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};
targetCells = {'freq'};
seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};
myCode = myCode(fileName,parametersCells,targetCells,seqOrder,4);
myCode.sizeOfCPG = 4;
% clear fileName parametersCells targetCells seqOrder

problemType = '4N_symm';

%% plot graph of perf over hidden neurons num
hiddenN = [4];
numOfRepeats = 20;
myCode.NN_Perf_over_HNnum(numOfRepeats,hiddenN,'train' );
pause(30);
myCode.NN_Perf_over_HNnum(numOfRepeats,hiddenN,'plot' );
myCode.NN_Perf_over_HNnum(numOfRepeats,hiddenN,'text' );
%% with NN:
myCode.disp_information = false;
myCode = myCode.Set('NN',[300],50);
myCode = myCode.trainNN(1,0);
myCode.plot_fit_data('NN',problemType);

% g = gpuDevice(1); % reset the GPU and clear memory
% reset(g); clear g

%% with the paper's MoE
myCode = myCode.Set('paper_MoE',10,5,0.005,0.995);
myCode = myCode.paper_MoE_train();
myCode.paper_MoE_plotPerf({'MSE_over_iter','regGraph','reg_with_color','reg_graph_from_NNtoolbox'})
myCode.plot_fit_data('paper_MoE',problemType);

%%
close all
myCode = myCode.Set('our_MoE',500,10,[100],1);
myCode.disp_information = true;
% myCode = myCode.my_MoE_train_collaboration();
myCode = myCode.my_MoE_train_Competetive('soft');
myCode.plot_fit_data('our_MoE',problemType);

%%
N=5;
training_epochs = 500;
myCode.disp_information = false;

mse_NN = zeros(N,1);
myCode = myCode.Set('NN',[4],training_epochs);
for i=1:N
    disp(['NN inter #',num2str(i)]);
    myCode = myCode.trainNN(0);
    mse_NN(i,1) = myCode.NN.MSE_test_perf;
    myCode = myCode.shuffle_samples('onlyTrainAndValid'); 
end

mse_MoE = zeros(N,1);
myCode = myCode.Set('our_MoE',training_epochs,2,[2],1);
for i=1:N
    disp(['MoE inter #',num2str(i)]);
%     myCode = myCode.my_MoE_train_Competetive('hard');
    myCode = myCode.my_MoE_train_collaboration();
    mse_MoE(i,1) = myCode.my_MoE_out.Moe_MSE_on_test;
    myCode = myCode.shuffle_samples('onlyTrainAndValid'); 
end

meanMse_NN = mean(mse_NN,1);
stdMse_NN = std(mse_NN,0,1);
meanMse_MoE = mean(mse_MoE,1);
stdMse_MoE = std(mse_MoE,0,1);

means = [meanMse_NN; meanMse_MoE];
stdevs = [stdMse_NN; stdMse_MoE];
Names = {' ',' ','   21 weights   ',' ',' '};
label_Y = 'perf_{MSE}';
graph_title = 'comparison between NN and MoE';
graph_legend = {'NN','MoE'};
myCode.plot_MSE_perf_with_stdev(means,stdevs,Names,label_Y,graph_title,graph_legend)
