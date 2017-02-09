
clear all; close all; clc

fileName = 'MatsRandomRes_2Neurons_change_only_a.mat';
parametersCells = {'tau','b','a','s'};
targetCells = {'freq'};
seqOrder = {'tau','b','s','_','a'}; %'_' - dont care parameter
myCode = myCode(fileName,parametersCells,targetCells,seqOrder);

myCode = myCode.Set('NN',[5,5]);
myCode = myCode.trainNN(1);