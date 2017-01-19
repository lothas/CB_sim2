clear all; clc;
genome_file = 'MatsuokaGenome.mat';
nAnkle = 1;%1; % Number of ankle torques
nHip = 0;   % Number of hip torques
maxAnkle = 20;   % Max ankle torque
maxHip = 20;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;
Mw = 10*ones(1,(2*N-1)*2*N);
mw = 0*Mw;
% %     % Final genome with tau_r + beta (constant tau_u/tau_v ratio) %%
    Keys = {'\tau_r', 'beta',        'amp',        '2neuron_general_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
                  1 ,      1,          2*N,                                2,        1 ,       2*N ,            0 };
    Range = {  0.02 ,    0.2,         mamp,                            [1,1],   -0.001 ,  -0.2*Mamp; % Min
               0.6  ,   10.0,         Mamp,                          [10,10],   0.001 ,   0.2*Mamp}; % Max

% % %     % 2neuron general specific range%%
% Mamp = [4,4];
% mamp = [4.01,4.01];
%     Keys = {'\tau_r', 'beta',        'amp',        '2neuron_general_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%                   1 ,      1,          2*N,                                2,        1 ,       2*N ,            0 };
%     Range = {  0.23 ,    2.4,         mamp,                            [1,1],   -0.001 ,  -0.2*Mamp; % Min
%                0.26  ,   2.6,         Mamp,                            [4,4],   0.001 ,   0.2*Mamp}; % Max
           
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
MML.tEnd = 15;

% % Set constant beta
% MML.Sim.Con.beta = 7;

nPlotSamples = 0; % 10;

% Turn off findpeaks warning
warning('off','signal:findpeaks:largeMinPeakHeight');


%% Phase 1 - Run lots of Matsuoka simulations with different parameters
filename1 = 'MatsRandomRes_2Neurons_general_16_01_2017_E_narrowRange.mat';
% filename1 = 'MatsRandomRes_test.mat';
nSamples = 1000;
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