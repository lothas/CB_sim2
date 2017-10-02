
%% NN testing Script:
% IMPORTANT: dont forget to load the right genome file and to uptade
% 'MatsuokaML.m' to the right rettings
% 

clear all; close all; clc;

%% Create genome (only if necessary)
genome_file = 'MatsuokaGenome_4Neuron_tagaLike.mat';
nAnkle = 1; % Number of ankle torques
nHip = 1;   % Number of hip torques
maxAnkle = 20;   % Max ankle torque
maxHip = 20; %8;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;
mw = 0*ones(1,4);

% CPG strucute: (ALSO Symm W_ij = W_ji)
%   H_F   H_E           % 
% 4 O-----O 3           %   
%    \    |             %   w = [0  , W12, 0  , 0  ; 
%     \   |             %        W21, 0  , w23, W24;
%      \  |             %        0  , 0  , 0  , W34;
%       \ |             %        0  , 0  , w43, 0  ;
% 1 O-----O 2           % w12=w21 = w1  
%  A_F    A_E           % w23 = w2
%                       % w24 = w3
%                       % w43=w34 = w4


Mw = 5*ones(1,4);
    % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
Keys = {'\tau_r', 'beta', 'amp_4n_symm',   '4neuron_taga_like', 'ks_\tau',     'ks_c_4n_symm', 'IC_matsuoka';
              1 ,      1,             1,                     4,        1 ,                 1 ,            0 };
Range = {  0.02 ,    0.2,             0,                    mw,      -10 ,               -0.1; % Min
           0.25 ,    2.5,        maxHip,                    Mw,       10 ,                0.1}; % Max

       
MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

clear all

%%

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau','b','c','w_{1}','w_{2}','w_{3}','w_{4}','k_tau','k_{c}'};

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;



