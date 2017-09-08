clear all; clc;
%%
genome_file = 'MatsuokaGenome_2Neuron_General.mat';
nAnkle = 1;%1; % Number of ankle torques
nHip = 0;   % Number of hip torques
maxAnkle = 5;%20;   % Max ankle torque
maxHip = 5;%20;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;


% %     % 2neuron general specific range%%
% % Large b Large W Ranges
% Mw = 10*ones(1,(2*N-1)*2*N);
% mw = 0*Mw;
% Keys = {'\tau_r', 'beta',        'amp',        '2neuron_general_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%               1 ,      1,            2,                                2,        1 ,       2 ,            0 };
% Range = {  0.02 ,      0,        [0,0],                               mw,   -0.001 ,  [-0.2,-0.2]; % Min
%            0.25 ,     10,        [5,5],                               Mw,    0.001 ,   [0.2,0.2]}; % Max
%    
% % Narrow b Large W Ranges
% Mw = 10*ones(1,(2*N-1)*2*N);
% mw = 0*Mw;
% Keys = {'\tau_r', 'beta',        'amp',        '2neuron_general_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%               1 ,      1,            2,                                2,        1 ,       2 ,            0 };
% Range = {  0.02 ,      0,        [0,0],                               mw,   -0.001 ,  [-0.2,-0.2]; % Min
%            0.25 ,     10,        [5,5],                               Mw,    0.001 ,   [0.2,0.2]}; % Max
% 
       % Narrow b Narrow W Ranges
Mw = 5*ones(1,(2*N-1)*2*N);
mw = 0*Mw;
Keys = {'\tau_r', 'beta',        'amp',        '2neuron_general_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
              1 ,      1,            2,                                2,        1 ,       2 ,            0 };
Range = {  0.02 ,      0,        [0,0],                               mw,   -0.001 ,  [-0.2,-0.2]; % Min
           0.25 ,     10,        [5,5],                               Mw,    0.001 ,   [0.2,0.2]}; % Max
    
       
MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

clear all
%% Initialize machine learning object for Matsuoka analysis
% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 2;

% % change tau_a/tau_r to 12 (instead of 5)
MML.Sim.Con.tau_ratio = 12;

% Turn off findpeaks warning
warning('off','signal:findpeaks:largeMinPeakHeight');


%% Phase 1 - Run lots of Matsuoka simulations with different parameters

% seqOrder = {'tau','b','c1','c2','w12','w21'};

N = 200000; % the number of samples
% CPG parameters:
tau_min = 0.02;     tau_max = 0.25;
tau = (tau_max-tau_min).*rand(1,N) + tau_min;

b_min = 0.2;     b_max = 2.5;
b = (b_max-b_min).*rand(1,N) + b_min;

c_min = 0;     c_max = 5;
c = (c_max-c_min).*rand(1,N) + c_min;

w12_min = 0;     w12_max = 5;
w12 = (w12_max-w12_min).*rand(1,N) + w12_min;

w21_min = 0;     w21_max = 5;
w21 = (w21_max-w21_min).*rand(1,N) + w21_min;

disp('start with the sim:');
parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
% for i=1:N
    disp(['at sim #',num2str(i)]);
    seq = [tau(1,i),b(1,i),c(1,i),c(1,i),w12(1,i),w21(1,i),0,0,0];
    [out, ~, signal] = MML.runSim(seq);
        % Prepare output:
    % Parameters
    results(i).seq = seq;
    results(i).tau = tau(1,i);
    results(i).b = b(1,i);
    results(i).c = c(1,i);
    results(i).w12 = w12(1,i);
    results(i).w21 = w21(1,i);
    results(i).x0 = out.x0;

    results(i).periods = out.periods;

    results(i).pos_work = out.pos_work;
    results(i).neg_work = out.neg_work;
    results(i).perError1 = out.perError1;
    results(i).perOK1 = out.perOK1;
    results(i).perError2 = out.perError2;
    results(i).perOK2 = out.perOK2;
    results(i).neuronActive = out.neuronActive;
    results(i).neuronOsc = out.neuronOsc;
    
end 
disp('sim end...');

header = sprintf('tau ratio is equal to 12 \n');
header = [header,sprintf('data is for 2N assymetric case \n')];
header = [header,sprintf('seq Order: \n')];
header = [header,sprintf('"tau","b","c_1","c_2","W_12","W_21" \n')];
header = [header,sprintf('b in range (0.2,2.5) \n')];
header = [header,sprintf('W_i in range (0,5) \n')];

save('MatsRandomRes_2Neurons_general_Narrow_b_Narrow_W_1.mat','results','header','MML');


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