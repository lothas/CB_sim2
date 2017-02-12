clear all; clc;
%%
genome_file = 'MatsuokaGenome.mat';
nAnkle = 1; % Number of ankle torques
nHip = 1;   % Number of hip torques
maxAnkle = 10;%20;   % Max ankle torque
maxHip = 10;%20;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;
Mw = 10*ones(1,6);
mw = 0*Mw;

% %     % 2neuron general specific range%%
Keys = {'\tau_r', 'beta',  'amp_4n_symm',           '4neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
              1 ,      1,              1,                                6,        1 ,       2*N ,            0 };
Range = {   0.1 ,      0,              0,                      0*ones(1,6),   -0.001 ,  -0.2*Mamp; % Min
            0.5 ,      5,       maxAnkle,                      5*ones(1,6),   0.001 ,   0.2*Mamp}; % Max
           
MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

%% Initialize machine learning object for Matsuoka analysis
clear all; clc;
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01; % 0.01
MML.tEnd = 30;

nPlotSamples = 0; % 10;

% Turn off findpeaks warning
warning('off','signal:findpeaks:largeMinPeakHeight');


%% Phase 1 - Run lots of Matsuoka simulations with different parameters
filename1 = 'MatsRandomRes_4Neurons_symm.mat';
% filename1 = 'MatsRandomRes_test.mat';
nSamples = 20000;
MML.runRandomSims(nSamples, filename1);

%% plotting
% load(filename1)
[out, ~, signal] = MML.runSim(results(randsample(1:nSamples,1)).seq);

figure;
subplot(2,1,1);
plot(signal.T,signal.X);
xlabel('time[sec]');    ylabel('X_i');
title('X_i over time');
subplot(2,1,2)
plot(signal.T,signal.signal(1,:));