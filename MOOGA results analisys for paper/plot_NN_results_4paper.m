%% NN results:
clear all; close all; clc

results_fileName = 'MatsRandomRes_4Neurons_4Paper_for_MOOGA_try.mat';

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;

NNs_4paper = NNs_4paper(results_fileName,MML);

%% plot histograms of  osc_param Vs n-osc_param:
paramName = 'tau';
NNs_4paper.plot_oscParam_vs_NoscParam_hist(paramName)