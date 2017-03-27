
clear all; close all; clc

%% 2N symmetric
% fileName = 'MatsRandomRes_2Neurons_change_only_a.mat';
% fileName = {'MatsRandomRes_2Neurons_symm_trainData_wide_range.mat',...
%     'MatsRandomRes_2Neurons_symm_testData.mat'};
fileName = {'MatsRandomRes_2Neurons_symm_trainData_narrow_range.mat',...
    'MatsRandomRes_2Neurons_symm_testData.mat'};
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
fileName = 'MatsRandomRes_31_01_2017.mat';
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
%% with NN:
myCode.disp_information = false;
myCode = myCode.Set('NN',[2],20);
myCode = myCode.trainNN(0);
myCode.plot_fit_data('NN',problemType);

%% with my_MoE
close all
competetiveflag = 3;
myCode = myCode.Set('our_MoE',2,2,[2],[2],5,competetiveflag);
myCode = myCode.my_MoE_train();
myCode.my_MoE_plot_train_perf();
% myCode.my_MoE_testNet(myCode.sampl_test,myCode.targ_test,myCode.my_MoE_out.expertsNN,...
%     myCode.my_MoE_out.gateNet,competetiveflag,1);
% myCode.plot_fit_data('my_MoE',problemType);

%% with the paper's MoE
myCode = myCode.Set('paper_MoE',10,5,0.005,0.995);
myCode = myCode.paper_MoE_train();
myCode.paper_MoE_plotPerf({'MSE_over_iter','regGraph','reg_with_color','reg_graph_from_NNtoolbox'})
myCode.plot_fit_data('paper_MoE',problemType);

%%
close all; clc

N = 15;
numOfExperts = 20;
numOfHiddenNeurons = 5;
numOfEpochs = 100;

myCode.disp_information = false;

train_err_NN = zeros(1,N);
valid_err_NN = zeros(1,N);
test_err_NN = zeros(1,N);

for i=1:N  
    myCode = myCode.shuffle_samples('completeShuffle');
    myCode = myCode.Set('NN',[numOfHiddenNeurons],numOfEpochs);
    myCode = myCode.trainNN(0);
    
    train_err_NN(1,i) = myCode.NN.net_perf.best_perf;
    valid_err_NN(1,i) = myCode.NN.net_perf.best_vperf;
    test_err_NN(1,i) = myCode.NN.net_perf.best_tperf;
end

competetiveflag = 1;
train_err_MoE_cmpFlag1 = zeros(1,N);
valid_err_MoE_cmpFlag1 = zeros(1,N);
test_err_MoE_cmpFlag1 = zeros(1,N);

for i=1:N
    myCode = myCode.shuffle_samples('completeShuffle');
    myCode = myCode.Set('our_MoE',numOfEpochs/5,numOfExperts,[numOfHiddenNeurons],[5],5,competetiveflag);
    myCode = myCode.my_MoE_train();

    [train_err_MoE_cmpFlag1(1,i),~] = myCode.NN_perf_calc(myCode.targ_train,myCode.my_MoE_out.out_from_train,0,0,'train');
    [valid_err_MoE_cmpFlag1(1,i),~] = myCode.NN_perf_calc(myCode.targ_valid,myCode.my_MoE_out.out_from_valid,0,0,'valid');
    [test_err_MoE_cmpFlag1(1,i),~] = myCode.NN_perf_calc(myCode.targ_test,myCode.my_MoE_out.out_from_test,0,0,'test');
end

competetiveflag = 3;
train_err_MoE_cmpFlag3 = zeros(1,N);
valid_err_MoE_cmpFlag3 = zeros(1,N);
test_err_MoE_cmpFlag3 = zeros(1,N);

for i=1:N
    myCode = myCode.shuffle_samples('completeShuffle');
    myCode = myCode.Set('our_MoE',numOfEpochs/5,numOfExperts,[numOfHiddenNeurons],[5],5,competetiveflag);
    myCode = myCode.my_MoE_train();

    [train_err_MoE_cmpFlag3(1,i),~] = myCode.NN_perf_calc(myCode.targ_train,myCode.my_MoE_out.out_from_train,0,0,'train');
    [valid_err_MoE_cmpFlag3(1,i),~] = myCode.NN_perf_calc(myCode.targ_valid,myCode.my_MoE_out.out_from_valid,0,0,'valid');
    [test_err_MoE_cmpFlag3(1,i),~] = myCode.NN_perf_calc(myCode.targ_test,myCode.my_MoE_out.out_from_test,0,0,'test');
end

means = [mean(train_err_NN,2),mean(valid_err_NN,2),mean(test_err_NN,2);
    mean(train_err_MoE_cmpFlag1,2),mean(valid_err_MoE_cmpFlag1,2),mean(test_err_MoE_cmpFlag1,2);
    mean(train_err_MoE_cmpFlag3,2),mean(valid_err_MoE_cmpFlag3,2),mean(test_err_MoE_cmpFlag3,2)];
stdevs = [std(train_err_NN,[],2),std(valid_err_NN,[],2),std(test_err_NN,[],2);
    std(train_err_MoE_cmpFlag1,[],2),std(valid_err_MoE_cmpFlag1,[],2),std(test_err_MoE_cmpFlag1,[],2);
    std(train_err_MoE_cmpFlag3,[],2),std(valid_err_MoE_cmpFlag3,[],2),std(test_err_MoE_cmpFlag3,[],2)];
Names = {' ','NN',' ','MoE "winner takes all"',' ','MoE more similar to the paper'};
label_Y = 'mean MSE perf';
graph_title = 'perf of NN Vs MoE';
graph_legend = {'train','validation','test'};
myCode.plot_MSE_perf_with_stdev(means',stdevs',Names,label_Y,graph_title,graph_legend);

disp(['mean NN MSE = ',num2str(means(1))]);
disp(['stdev NN MSE = ',num2str(stdevs(1))]);

disp(['mean MoE MSE 1st method = ',num2str(means(2))]);
disp(['stdev MoE MSE 1st method = ',num2str(stdevs(2))]);

disp(['mean MoE MSE 3rd method = ',num2str(means(3))]);
disp(['stdev MoE MSE 3rd method = ',num2str(stdevs(3))]);


%% without the 1st method plot
close all; clc

N = 20;
numOfExperts = 100;
numOfHiddenNeurons = 2;
numOfEpochs = 100;

myCode.disp_information = false;

train_err_NN = zeros(1,N);
valid_err_NN = zeros(1,N);
test_err_NN = zeros(1,N);

for i=1:N  
    myCode = myCode.shuffle_samples('completeShuffle');
    myCode = myCode.Set('NN',[numOfHiddenNeurons],numOfEpochs);
    myCode = myCode.trainNN(0);
    
    train_err_NN(1,i) = myCode.NN.net_perf.best_perf;
    valid_err_NN(1,i) = myCode.NN.net_perf.best_vperf;
    test_err_NN(1,i) = myCode.NN.net_perf.best_tperf;
end

competetiveflag = 3;
train_err_MoE_cmpFlag3 = zeros(1,N);
valid_err_MoE_cmpFlag3 = zeros(1,N);
test_err_MoE_cmpFlag3 = zeros(1,N);

for i=1:N
    myCode = myCode.shuffle_samples('completeShuffle');
    myCode = myCode.Set('our_MoE',numOfEpochs/5,numOfExperts,[numOfHiddenNeurons],[5],5,competetiveflag);
    myCode = myCode.my_MoE_train();

    [train_err_MoE_cmpFlag3(1,i),~] = myCode.NN_perf_calc(myCode.targ_train,myCode.my_MoE_out.out_from_train,0,0,'train');
    [valid_err_MoE_cmpFlag3(1,i),~] = myCode.NN_perf_calc(myCode.targ_valid,myCode.my_MoE_out.out_from_valid,0,0,'valid');
    [test_err_MoE_cmpFlag3(1,i),~] = myCode.NN_perf_calc(myCode.targ_test,myCode.my_MoE_out.out_from_test,0,0,'test');
end

means = [mean(train_err_NN,2),mean(valid_err_NN,2),mean(test_err_NN,2);
    mean(train_err_MoE_cmpFlag3,2),mean(valid_err_MoE_cmpFlag3,2),mean(test_err_MoE_cmpFlag3,2)];
stdevs = [std(train_err_NN,[],2),std(valid_err_NN,[],2),std(test_err_NN,[],2);
    std(train_err_MoE_cmpFlag3,[],2),std(valid_err_MoE_cmpFlag3,[],2),std(test_err_MoE_cmpFlag3,[],2)];
Names = {' ',' ','NN',' ',' ',' ',' ',' ','MoE more similar to the paper',' ',' '};
label_Y = 'mean MSE perf';
graph_title = 'perf of NN Vs MoE';
graph_legend = {'train','validation','test'};
myCode.plot_MSE_perf_with_stdev(means',stdevs',Names,label_Y,graph_title,graph_legend);

disp(['mean NN MSE = ',num2str(means(1,3))]);
disp(['stdev NN MSE = ',num2str(stdevs(1,3))]);

disp(['mean MoE MSE 3rd method = ',num2str(means(2,3))]);
disp(['stdev MoE MSE 3rd method = ',num2str(stdevs(2,3))]);

%% comparing to various #experts MoE
close all; clc

N = 10;
numOfExperts = [5,10,20,50];
numOfHiddenNeurons = 5;
numOfEpochs = 50;

myCode.disp_information = false;

train_err_NN = zeros(1,N);
valid_err_NN = zeros(1,N);
test_err_NN = zeros(1,N);

for i=1:N  
    myCode = myCode.shuffle_samples('completeShuffle');
    myCode = myCode.Set('NN',[numOfHiddenNeurons],numOfEpochs);
    myCode = myCode.trainNN(0);
    
    train_err_NN(1,i) = myCode.NN.net_perf.best_perf;
    valid_err_NN(1,i) = myCode.NN.net_perf.best_vperf;
    test_err_NN(1,i) = myCode.NN.net_perf.best_tperf;
end

mean_MSE_NN = [mean(train_err_NN,2),mean(valid_err_NN,2),mean(test_err_NN,2)];
stdev_MSE_NN = [std(train_err_NN,[],2),std(valid_err_NN,[],2),std(test_err_NN,[],2)];

competetiveflag = 3;
numOfSims = length(numOfExperts);
train_err_MoE_cmpFlag3 = zeros(numOfSims,N);
valid_err_MoE_cmpFlag3 = zeros(numOfSims,N);
test_err_MoE_cmpFlag3 = zeros(numOfSims,N);

for j=1:numOfSims
    for i=1:N
        myCode = myCode.shuffle_samples('completeShuffle');
        myCode = myCode.Set('our_MoE',numOfEpochs/5,numOfExperts(1,j),[numOfHiddenNeurons],[5],5,competetiveflag);
        myCode = myCode.my_MoE_train();

        [train_err_MoE_cmpFlag3(j,i),~] = myCode.NN_perf_calc(myCode.targ_train,myCode.my_MoE_out.out_from_train,0,0,'train');
        [valid_err_MoE_cmpFlag3(j,i),~] = myCode.NN_perf_calc(myCode.targ_valid,myCode.my_MoE_out.out_from_valid,0,0,'valid');
        [test_err_MoE_cmpFlag3(j,i),~] = myCode.NN_perf_calc(myCode.targ_test,myCode.my_MoE_out.out_from_test,0,0,'test');
    end
end
mean_MSE_MoE = [mean(train_err_MoE_cmpFlag3,2),...
    mean(valid_err_MoE_cmpFlag3,2),...
    mean(test_err_MoE_cmpFlag3,2)];
stdev_MSE_MoE = [std(train_err_MoE_cmpFlag3,[],2),...
    std(valid_err_MoE_cmpFlag3,[],2),...
    std(test_err_MoE_cmpFlag3,[],2)];

means = [mean_MSE_NN;
    mean_MSE_MoE];
stdevs = [stdev_MSE_NN;
    stdev_MSE_MoE];
Names = {' ','NN',' ','MoE 2exp',' ','MoE 10exp',' ','MoE 20exp',' ','MoE 50exp',' '};
label_Y = 'mean MSE perf';
graph_title = 'perf of NN Vs MoE';
graph_legend = {'train','validation','test'};
myCode.plot_MSE_perf_with_stdev(means',stdevs',Names,label_Y,graph_title,graph_legend);

disp(['mean NN MSE = ',num2str(means(1,3))]);
disp(['stdev NN MSE = ',num2str(stdevs(1,3))]);

for j=1:numOfSims
    disp(['mean MoE MSE with ',num2str(numOfExperts(1,j)),'experts = ',num2str(means(j+1,3))]);
    disp(['stdev MoE MSE with ',num2str(numOfExperts(1,j)),'experts = ',num2str(stdevs(j+1,3))]);
end




