
clear all; close all; clc

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

%%

hiddenSize = 4;
autoenc = trainAutoencoder(myCode.sampl_train,hiddenSize,...
    'DecoderTransferFunction','purelin');
% plotWeights(autoenc)

% autoenc = trainAutoencoder(myCode.sampl_train,hiddenSize,...
%     'L2WeightRegularization',0.001,...
%     'SparsityRegularization',4,...
%     'SparsityProportion',0.05,...
%     'DecoderTransferFunction','purelin');
net = network(autoenc);

inputPred = net(myCode.sampl_train);

[C,S,B] = regression(myCode.sampl_train,inputPred,'one');

plotregression(myCode.sampl_train,inputPred);

%%
hiddenSize = 6;
autoenc = trainAutoencoder(myCode.sampl_train,hiddenSize,...
    'DecoderTransferFunction','purelin');
features1 = encode(autoenc,myCode.sampl_train);
hiddenSize = 4;
autoenc2 = trainAutoencoder(features1,hiddenSize,...
    'DecoderTransferFunction','purelin',...
    'ScaleData',false);
deepnet = stack(autoenc,autoenc2);

inputPred = deepnet(myCode.sampl_train);

inputs = myCode.sampl_test;
targets = myCode.targ_test;
net2 = feedforwardnet(5);
net2 = train(net2,inputs,targets);

[C,S,B] = regression(myCode.sampl_train,inputPred,'one');

plotregression(myCode.sampl_train,inputPred);
%%
inputPredTest = net(myCode.sampl_test);

[C,S,B] = regression(myCode.sampl_test,inputPredTest,'one');

plotregression(myCode.sampl_test,inputPredTest);

%%
inputs = myCode.sampl_test;
targets = myCode.targ_test;
features1 = encode(autoenc,inputs);

net2 = feedforwardnet(10);
net2 = train(net2,features1,targets);

net3 = feedforwardnet(10);
net3 = train(net3,myCode.sampl_train,myCode.targ_train);
plotregression(myCode.targ_test,net3(myCode.sampl_test),'Test');