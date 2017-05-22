
% resume work on NN for the paper script

close all; clear all; clc;

%% NN input outputs: (option 3)
parametersCells = {'freq_desired','b',...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};
targetCells = {'tau'};

%% NN input outputs: (option 5)
parametersCells = {'freq_desired',...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};
targetCells = {'b'};

%% NN input outputs: (option 9)
parametersCells = {'freq_desired',...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};
targetCells = {'tau','b'};

% note: not sure is would run with two outputs

%% File to load
% fileName = 'MatsRandomRes_test.mat';
fileName = 'MatsRandomRes_2-4_02_2017.mat';
seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};
problemType = '4N_symm';

%% Load data 
myCode = myCode(fileName,parametersCells,targetCells,seqOrder,4);
myCode.sizeOfCPG = 4;

%% train with NN: (general training)
myCode.disp_information = false;
myCode = myCode.Set('NN',[2],50);
myCode = myCode.trainNN(1,0);
myCode.plot_fit_data('NN',problemType);
% g = gpuDevice(1); % reset the GPU and clear memory
% reset(g); clear g

%% train with NN: (graph: perf over hidNum)
HiddenN = [8,10,12,14,16,18,20,50,100,300];
% HiddenN = [3,5];
NumOfRepeats = 5;
Smpl2HidN_ratio_const = true;
myCode.NN_Perf_over_HNnum(NumOfRepeats,HiddenN,'train',Smpl2HidN_ratio_const);
myCode.NN_Perf_over_HNnum(NumOfRepeats,HiddenN,'plot',Smpl2HidN_ratio_const);

%% train with NN: (graph: perf over num of samples)
NumOfRepeats = 10;
HiddenN = 10;
dataPointsNum = 1000 * [15,20,25,30,40,50,60,120,150,200,225];
% % dataPointsNum = 1000 * [15,20];
NNinput = horzcat(myCode.sampl_train,myCode.sampl_valid);
NNtarg = horzcat(myCode.targ_train,myCode.targ_valid);

myCode.NN_Perf_over_sampl_num(NumOfRepeats,HiddenN,dataPointsNum,...
    NNinput,NNtarg,'train')

myCode.NN_Perf_over_sampl_num(NumOfRepeats,HiddenN,dataPointsNum,...
    NNinput,NNtarg,'plot')

%% TODO:
% add code to save data as .cvs for tensorflow use

%% TODO2:
% use 500000 samples for training