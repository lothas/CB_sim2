
%% Initialize machine learning object for Matsuoka analysis
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01; % 0.01
MML.tEnd = 15; % 15

% % Set constant beta
% MML.Sim.Con.beta = 7;

nPlotSamples = 0; % 10;

% Turn off findpeaks warning
warning('off','signal:findpeaks:largeMinPeakHeight');


%% Phase 1 - Run lots of Matsuoka simulations with different parameters
filename1 = 'MatsRandomRes_2Neurons_10_01_2017_O.mat';
% filename1 = 'MatsRandomRes_test.mat';
nSamples = 70000;
MML.runRandomSims(nSamples, filename1);
