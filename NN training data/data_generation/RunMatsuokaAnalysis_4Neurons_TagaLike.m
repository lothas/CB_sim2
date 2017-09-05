
close all; clc; clear all;

%%
genome_file = 'MatsuokaGenome_4Neuron_tagaLike.mat';
nAnkle = 1; % Number of ankle torques
nHip = 1;   % Number of hip torques
maxAnkle = 20;   % Max ankle torque
maxHip = 8;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;
mw = 0*ones(1,4);

% CPG strucute: (ALSO Symm W_ij = W_ji)
%   H_F   H_E           % 
% 4 O-----O 3           %   
%    \    |             %   w = [0  , W12, 0  , 0  ; 
%     \   |             %        W21, 0  , w23, W24;
%      \  |             %        0  , w32, 0  , W34;
%       \ |             %        0  , W42, w43, 0  ;
% 1 O-----O 2           % w12=w21 = w1  
%  A_F    A_E           % w23=w32 = w2
%                       % w24=w42 = w3
%                       % w43=w34 = w4

%     %%%%%%%%%%%% For the 4-neuron case!!!
% % large b, large W
% Mw = 10*ones(1,4);
% %     % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
%     Keys = {'\tau_r', 'beta', 'amp',   '4neuron_taga_like', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%                   1 ,      1,    4 ,           4,        1 ,         4 ,            0 };
%     Range = {  0.02 ,    0.2,  mamp,          mw,      -10 ,  -0.1*Mamp; % Min
%                0.25 ,   10.0,  Mamp,          Mw,       10 ,   0.1*Mamp}; % Max
%          

% %%%%%%%%%%% For the 4-neuron case!!!
% % Narrow b, Narrow W
% Mw = 5*ones(1,4);
%     % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
% Keys = {'\tau_r', 'beta', 'amp',   '4neuron_taga_like', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%               1 ,      1,    4 ,           4,        1 ,         4 ,            0 };
% Range = {  0.02 ,    0.2,  mamp,          mw,      -10 ,  -0.1*Mamp; % Min
%            0.25 ,    2.5,  Mamp,          Mw,       10 ,   0.1*Mamp}; % Max

%%%%%%%%%%%% For the 4-neuron case!!!
% Narrow b, Large W
Mw = 10*ones(1,4);
%     % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
Keys = {'\tau_r', 'beta', 'amp',   '4neuron_taga_like', 'ks_\tau',     'ks_c', 'IC_matsuoka';
              1 ,      1,    4 ,           4,        1 ,         4 ,            0 };
Range = {  0.02 ,    0.2,  mamp,          mw,      -10 ,  -0.1*Mamp; % Min
           0.25 ,    2.5,  Mamp,          Mw,       10 ,   0.1*Mamp}; % Max

MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

clear all
%%
% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;

% % change tau_a/tau_r to 12 (instead of 5)
MML.Sim.Con.tau_ratio = 12;
% MML.Gen.Range(2,2) = 2.5; % the class will filter genes that are not in the new range.
%% Train data:

N = 400000; % the number of samples
% % % % CPG parameters:

tau_min = 0.02;     tau_max = 0.25;
tau = (tau_max-tau_min).*rand(1,N) + tau_min;

% b_min = 0.2;     b_max = 10;
b_min = 0.2;     b_max = 2.5;
b = (b_max-b_min).*rand(1,N) + b_min;

c_hip_min = 0;     c_hip_max = 8;
c_hip = (c_hip_max-c_hip_min).*rand(2,N) + c_hip_min;
c_ankle_min = 0;     c_ankle_max = 20;
c_ankle = (c_ankle_max-c_ankle_min).*rand(2,N) + c_ankle_min;
c = [c_ankle;c_hip];

W_min = 0;     W_max = 10;
% W_min = 0;     W_max = 5;
W = (W_max-W_min).*rand(4,N) + W_min;

ks_tau_min = -10;     ks_tau_max = 10;
ks_tau = (ks_tau_max-ks_tau_min).*rand(1,N) + ks_tau_min;

ks_c_hip_min = -0.8;     ks_c_hip_max = 0.8;
ks_c_hip = (ks_c_hip_max-ks_c_hip_min).*rand(2,N) + ks_c_hip_min;
ks_c_ankle_min = -2;     ks_c_ankle_max = 2;
ks_c_ankle = (ks_c_ankle_max-ks_c_ankle_min).*rand(2,N) + ks_c_ankle_min;
ks_c = [ks_c_ankle;ks_c_hip];

tic
t_cur = tic;

disp('start with the sim:');
parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
% for i=1:N
    disp(['at sim #',num2str(i)]);
    seq = [tau(1,i),b(1,i),c(:,i)',W(:,i)',ks_tau(1,i),ks_c(:,i)'];
    [out, sim, signal] = MML.runSim(seq);
        % Prepare output:
    % Parameters
    results(i).seq = seq;
    results(i).b = sim.Con.beta;
    results(i).c = sim.Con.Amp0;
    results(i).Worig = sim.Con.wex;
    results(i).W = sim.Con.W;
    results(i).Tr = sim.Con.tau;
    results(i).Ta = sim.Con.tav;
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

t_elapsed = toc(t_cur);
avg_sim_time = t_elapsed/N;
disp(['avg sim time is ',num2str(avg_sim_time),' [sec]']);

header = sprintf('tau ratio is equal to 12 \n');
header = [header,sprintf('data is for 4N TagaLike case \n')];
header = [header,sprintf('seq Order: \n')];
header = [header,sprintf('"tau","b","c_1","c_2","c_3","c_4" \n')];
header = [header,sprintf('"w_{1}","w_{2}","w_{3}","w_{4}" \n')];
header = [header,sprintf('b in range (0.2,2.5) \n')];
header = [header,sprintf('W_i in range (0,10) \n')];

save('MatsRandomRes_4Neurons_TagaLike_Narrow_b_Large_W_2.mat','results','MML','header');

%% plot example:
clc; close all

N = length(results);
rand_id = randsample(1:N,1);

[out, ~, signal] = MML.runSim(results(rand_id).seq);

figure;
subplot(2,1,1);
plot(signal.T,signal.X);
xlabel('time[sec]');    ylabel('X_i');
title({'X_i over time',...
    ['id #',num2str(rand_id),...
    '    periods: ',...
    num2str(results(rand_id).periods(1)),...
    ' ',num2str(results(rand_id).periods(2))]});
subplot(2,1,2)
plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');

clear out signal