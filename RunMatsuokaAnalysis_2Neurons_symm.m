
close all; clc; clear all;

%%
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
% %     % 2neuron symmetric specific range%%
Keys = {'\tau_r', 'beta',     'amp_2n',        '2neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
              1 ,      1,          2*N,                             1,        1 ,       2*N ,            0 };
Range = {  0.02 ,    1.1,         mamp,                             1,   -0.001 ,  -0.2*Mamp; % Min
           0.6  ,    8.0,         Mamp,                             6,   0.001 ,   0.2*Mamp}; % Max

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
MML.tEnd = 50; % 15
MML.nNeurons = 2;
%% Train data:
N = 100000; % the number of samples
% CPG parameters:
% tau_min = 0.4;     tau_max = 0.6;
tau_min = 0.02;     tau_max = 0.6;
tau = (tau_max-tau_min).*rand(1,N) + tau_min;

% b_min = 2;     b_max = 3;
b_min = 1.1;     b_max = 5;
b = (b_max-b_min).*rand(1,N) + b_min;

c_min = 2.5;     c_max = 7;
c = (c_max-c_min).*rand(1,N) + c_min;

a_min = 1;     a_max = 6;
a = (a_max-a_min).*rand(1,N) + a_min;

disp('start with the sim:');
parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
% for i=1:N
    disp(['at sim #',num2str(i)]);
    seq = [tau(1,i),b(1,i),c(1,i),0,a(1,i),0,0,0];
    [out, ~, signal] = MML.runSim(seq);
        % Prepare output:
    % Parameters
    results(i).seq = seq;
    results(i).tau = tau(1,i);
    results(i).b = b(1,i);
    results(i).c = c(1,i);
    results(i).a = a(1,i);
    results(i).x0 = out.x0;

    % Results- caculate perdiods using different methods:
    results(i).periods = out.periods;
    % results(i).periods_Rea = out.periods_Rea;
    
    t = signal.T;
    x = signal.signal(1,:);
    
    [periodsFFT,sine_coef,cos_coef,a0,~] = ...
        MML.processResults_LSQ(x,t,0);
    results(i).periods_LSQ = periodsFFT;
    results(i).bias_coef = a0;
    results(i).sine_coef = sine_coef;
    results(i).cos_coef = cos_coef;
    
    results(i).pos_work = out.pos_work;
    results(i).neg_work = out.neg_work;
    results(i).perError1 = out.perError1;
    results(i).perOK1 = out.perOK1;
    results(i).perError2 = out.perError2;
    results(i).perOK2 = out.perOK2;
    results(i).neuronActive = out.neuronActive;
    results(i).neuronOsc = out.neuronOsc;
    
    T = 5*tau(1,i);
    wn = (1/T)*sqrt(((tau(1,i)+T)*b(1,i)-tau(1,i)*a(1,i))/(tau(1,i)*a(1,i)));
    results(i).period_matsuoka_est = 2*pi./wn;
    
    Kn = (tau(1,i)+T)/(T*a(1,i));
    results(i).amp_matsuoka_est = c(1,i) / (2*Kn-1+(2/pi)*(a(1,i)+b(1,i))*asin(Kn))
%     clear t x periodsFFT sine_coef cos_coef a0
end 
disp('sim end...');

save('MatsRandomRes_2Neurons_symm_test_for_FFT_new_3.mat','results');


%% Test data:
% the test data has only change of 'a'
N = 1000; % the number of samples
% CPG parameters:
tau_min = 0.02;     tau_max = 0.6;
tau = 0.5.*ones(1,N);%(tau_max-tau_min).*rand(1,N) + tau_min;

b_min = 1.1;     b_max = 5;
b = 2.5.*ones(1,N);%(b_max-b_min).*rand(1,N) + b_min;

c_min = 0;     c_max = 10;
c = 5.*ones(1,N);

a_min = 1;     a_max = 6;
a = (a_max-a_min).*rand(1,N) + a_min;

disp('start with the sim:');
parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
    disp(['at sim #',num2str(i)]);
    seq = [tau(1,i),b(1,i),c(1,i),0,a(1,i),0,0,0];
    [out, ~, ~] = MML.runSim(seq);
        % Prepare output:
    % Parameters
    results(i).seq = seq;
    results(i).tau = tau(1,i);
    results(i).b = b(1,i);
    results(i).c = c(1,i);
    results(i).a = a(1,i);
    results(i).x0 = out.x0;

    % Results
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

save('MatsRandomRes_2Neurons_symm_testData2.mat','results');