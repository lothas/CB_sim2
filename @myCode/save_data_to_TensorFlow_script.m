% resume work on NN for the paper script TENSORFLOW

close all; clear all; clc;

%% NN input outputs: (option 3)
parametersCells = {'freq','b',...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};
targetCells = {'tau'};

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

%% Get input and output data:
train_and_valid_sampl = horzcat(myCode.sampl_train,myCode.sampl_valid);
train_and_valid_targ = horzcat(myCode.targ_train,myCode.targ_valid);
test_sampl = myCode.sampl_test;
test_targ = myCode.targ_test;

train_and_valid = vertcat(train_and_valid_sampl,train_and_valid_targ);
test = vertcat(test_sampl,test_targ);

%% make files names:
dataFileName_trval = ['dataForTensorFlow_4N_trainValid_',datestr(now,'mm_dd_hh_MM'),'.csv'];
dataFileName_test = ['dataForTensorFlow_4N_test',datestr(now,'mm_dd_hh_MM'),'.csv'];
headerFileName = ['dataForTensorFlow_4N_header',datestr(now,'mm_dd_hh_MM'),'.txt'];
%% save data file as CVS
csvwrite(dataFileName_trval,train_and_valid);
csvwrite(dataFileName_test,test);

%% make header file:
source_of_data = ['data file is: ',fileName,'.mat'];
% NOTE: make sure the you change this order according to the relevant
%       sequence!
seq_order = ['the order of the seq is: \n','1)freq \n','2)b \n',...
    '3) w_{12} \n','4) w_{13} \n','5) w_{14} \n',...
    '6) w_{21} \n','7) w_{23} \n','8) w_{24} \n',...
    '9) w_{31} \n','10) w_{32} \n','11) w_{34} \n',...
    '12) w_{41} \n','13) w_{42} \n','14) w_{43} \n','15) tau'];

fid = fopen(headerFileName,'wt');
fprintf(fid,[source_of_data,'\n']);
fprintf(fid,[seq_order,'\n']);
fclose(fid);



