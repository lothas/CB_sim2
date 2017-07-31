
close all; clc; clear all;

%%
genome_file = 'MatsuokaGenome.mat';
nAnkle = 1; % Number of ankle torques
nHip = 1;   % Number of hip torques
maxAnkle = 20;   % Max ankle torque
maxHip = 8;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;
Mw = 10*ones(1,12);
mw = 0*Mw;

    %%%%%%%%%%%% For the 4-neuron case!!!
%     % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
    Keys = {'\tau_r', 'beta', 'amp',   'weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
                  1 ,      1,    4 ,          12,        1 ,         4 ,            0 };
    Range = {  0.02 ,    0.2,  mamp,          mw,      -10 ,  -0.1*Mamp; % Min
               0.25 ,   10.0,  Mamp,          Mw,       10 ,   0.1*Mamp}; % Max
         
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
%% Train data:

N = 200000; % the number of samples
% % % % CPG parameters:
% seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
%     'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
%     'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};

tau_min = 0.02;     tau_max = 0.25;
tau = (tau_max-tau_min).*rand(1,N) + tau_min;

b_min = 0.2;     b_max = 10;
b = (b_max-b_min).*rand(1,N) + b_min;

c_hip_min = 0;     c_hip_max = 8;
c_hip = (c_hip_max-c_hip_min).*rand(2,N) + c_hip_min;
c_ankle_min = 0;     c_ankle_max = 20;
c_ankle = (c_ankle_max-c_ankle_min).*rand(2,N) + c_ankle_min;
c = [c_ankle;c_hip];

W_min = 0;     W_max = 10;
W = (W_max-W_min).*rand(12,N) + W_min;

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

    % Results- caculate perdiods using different methods:
    results(i).periods = out.periods;
    
%     [periods_LSQ1,~,~,~,~] = ...
%         MML.processResults_LSQ(signal.signal(1,:),signal.T,0);
%     [periods_LSQ2,~,~,~,~] = ...
%         MML.processResults_LSQ(signal.signal(2,:),signal.T,0);
%     results(i).periods_LSQ = [periods_LSQ1(1,1);periods_LSQ2(1,1)];
%     [periodsFFT,sine_coef,cos_coef,a0,~] = ...
%         MML.processResults_LSQ(x,t,0);
%     results(i).bias_coef = a0;
%     results(i).sine_coef = sine_coef;
%     results(i).cos_coef = cos_coef;
    
    results(i).pos_work = out.pos_work;
    results(i).neg_work = out.neg_work;
    results(i).perError1 = out.perError1;
    results(i).perOK1 = out.perOK1;
    results(i).perError2 = out.perError2;
    results(i).perOK2 = out.perOK2;
    results(i).neuronActive = out.neuronActive;
    results(i).neuronOsc = out.neuronOsc;
    
    % check Matsuoka conditions:
%     results(i).cond0 = checkCond_0(results(i).seq,seqOrder);
%     
%     results(i).cond1 = checkCond_1(results(i).seq,seqOrder);
%     
%     results(i).cond2 = checkCond_2(results(i).seq,seqOrder,cond1);
    
%     disp(['at i=',num2str(i),...
%         '  cond0=',num2str(results(i).cond0),...
%         '  cond1=[',num2str(results(i).cond1),']',...
%         '  cond2=',num2str(results(i).cond2)]);

end 
disp('sim end...');

t_elapsed = toc(t_cur);
avg_sim_time = t_elapsed/N;
disp(['avg sim time is ',num2str(avg_sim_time),' [sec]']);

save('MatsRandomRes_4Neurons_4Paper_7.mat','results');

%% Phase 2 - Re-run simulations that converged outside the desired range,
% this time with scaled temporal parameters
filename2 = 'MatsScaledRes.mat';
load(filename1,'results');
periods = horzcat(results(:).periods);
reDo_ids = zeros(1, length(results));
reDo_ids(~isnan(periods)) = 1;
reDo_ids(periods >= MML.perLimOut(1) &...
    periods <= MML.perLimOut(2)) = 0;

inputData = results(logical(reDo_ids));
inputPeriods = periods(logical(reDo_ids));
MML.runScaledSims(inputData, inputPeriods, filename2);

%% plot example:
% load(filename1)
[out, ~, signal] = MML.runSim(results(randsample(1:nSamples,1)).seq);

figure;
subplot(2,1,1);
plot(signal.T,signal.X);
xlabel('time[sec]');    ylabel('X_i');
title('X_i over time');
subplot(2,1,2)
plot(signal.T,signal.signal(1,:));
