
close all; clear all; clc

%%
genome_file = 'MatsuokaGenome_2Neuron_Symm.mat';
nAnkle = 1;%1; % Number of ankle torques
nHip = 1;   % Number of hip torques
maxAnkle = 10;   % Max ankle torque
maxHip = 10;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;
% %     % 2neuron symmetric specific range%%

       % % % Narrow b Narrow W Narrow tau
Mw = 5*ones(1,(2*N-1)*2*N);
mw = 0*Mw;
Keys = {'\tau_r', 'beta','amp_2n_same_inputs',    '2neuron_symm_weights', 'ks_\tau',     'ks_c_2n_symm', 'IC_matsuoka';
              1 ,      1,                   2,                         1,        1 ,          1,            0 };
Range = {  0.02 ,    0.2,               [0,0],                         0,   -0.001 ,       -0.2; % Min
           0.10  ,     5,             [10,10],                         5,    0.001 ,       0.2}; % Max

       
MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

clear all

%%
% GA_try_2N_Symm_Matsuoka('GA + NN_classi',[]);

GA_try_2N_Symm_Matsuoka('GA + NN_classi + rescale',[]);