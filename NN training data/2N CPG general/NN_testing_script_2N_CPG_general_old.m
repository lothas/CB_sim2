
%% NN testing Script:
% IMPORTANT: dont forget to load the right genome file and to uptade
% 'MatsuokaML.m' to the right rettings
% 

clear all; close all; clc;

%% Create genome (only if necessary)
genome_file = 'MatsuokaGenome_2Neuron_General.mat';
nAnkle = 1;%1; % Number of ankle torques
nHip = 0;   % Number of hip torques
maxAnkle = 5;%20;   % Max ankle torque
maxHip = 5;%20;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;


% %     % 2neuron general specific range%%
% Large b Large W Ranges
Mw = 10*ones(1,(2*N-1)*2*N);
mw = 0*Mw;
Keys = {'\tau_r', 'beta',        'amp',        '2neuron_general_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
              1 ,      1,            2,                                2,        1 ,       2 ,            0 };
Range = {  0.02 ,      0,        [0,0],                               mw,   -0.001 ,  [-0.2,-0.2]; % Min
           0.25 ,     10,        [5,5],                               Mw,    0.001 ,   [0.2,0.2]}; % Max
   
% % Narrow b Large W Ranges
% Mw = 10*ones(1,(2*N-1)*2*N);
% mw = 0*Mw;
% Keys = {'\tau_r', 'beta',        'amp',        '2neuron_general_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%               1 ,      1,            2,                                2,        1 ,       2 ,            0 };
% Range = {  0.02 ,      0,        [0,0],                               mw,   -0.001 ,  [-0.2,-0.2]; % Min
%            0.25 ,     10,        [5,5],                               Mw,    0.001 ,   [0.2,0.2]}; % Max

%        % Narrow b Narrow W Ranges
% Mw = 5*ones(1,(2*N-1)*2*N);
% mw = 0*Mw;
% Keys = {'\tau_r', 'beta',        'amp',        '2neuron_general_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%               1 ,      1,            2,                                2,        1 ,       2 ,            0 };
% Range = {  0.02 ,      0,        [0,0],                               mw,   -0.001 ,  [-0.2,-0.2]; % Min
%            0.25 ,     10,        [5,5],                               Mw,    0.001 ,   [0.2,0.2]}; % Max
%     
       
MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

clear all
%%

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau' ,'b', 'c_1', 'c_2', 'W_12','W_21'};
% "NR" - not relevnt param 

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 2;

% % change tau_a/tau_r to 12 (instead of 5)
MML.Sim.Con.tau_ratio = 12;

% file name for uploading:L
results_fileName = {'MatsRandomRes_2Neurons_general_Large_b_Large_W.mat'};
% results_fileName = {'MatsRandomRes_2Neurons_general_Narrow_b_Narrow_W_1.mat'};
% results_fileName = {'MatsRandomRes_2Neurons_general_Narrow_b_Large_W_1.mat'};
%% Load data:
load(results_fileName{1,1},'results','header');
disp('data file information:');
disp(header);
% if I have more than 1 data file:
for i=2:numel(results_fileName)
    data = load(results_fileName{1,i},'results');
    results = [results, data.results]; %#ok<AGROW>
end

clear data i

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
    num2str(results(rand_id).periods(1))]});
subplot(2,1,2)
plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');

clear out signal N rand_id

%% get and filter periods:

% get CPG periods:
periods = horzcat(results(:).periods);

% Filter CPG's where not both signals oscillating:
osc_ids_temp = ~isnan(periods);
osc_ids_temp = osc_ids_temp(1,:);
disp(['Number of non-osc CPGs: ',num2str(sum(~osc_ids_temp))]);

% % check that all of the parameters are in the genome range:
seq = (vertcat(results(:).seq))';
ids_in_genome_range = true(1,size(seq,2));
for n=1:MML.Gen.Length
    ids_temp = (seq(n,:) > MML.Gen.Range(1,n)) &...
        (seq(n,:) < MML.Gen.Range(2,n));
    ids_in_genome_range = ids_in_genome_range & ids_temp;
end
disp(['Number of CPGs with parameters not in range: ',...
    num2str(sum(~ids_in_genome_range))]);

num_of_osc_ids_exluded_param_range = sum(osc_ids_temp);
osc_ids = osc_ids_temp & ids_in_genome_range;

osc_inRange_ids = osc_ids &...
    ( (periods(1,:) > MML.perLimOut(1,1)) &...
    (periods(1,:) < MML.perLimOut(1,2)) );

num_of_osc_ids = sum(osc_ids);
num_of_inRange_ids = sum(osc_inRange_ids);

%% keep the good seq and periods
seq = seq(:,osc_ids);
periods = periods(:,osc_ids);

%% Prepare NN inputs and outputs:
% input_names = {'b','tau','W_12','W_21'};
% output_names = {'periods'};

input_names = {'periods','tau','W_12','W_21'};
output_names = {'b'};


[sampl,targ] = ...
    prepare_NN_data(input_names,output_names,...
    seqOrder,seq,periods);

%% Neural Network:
architecture = [20];

net = fitnet(architecture);
net.trainFcn = 'trainbr';
net.trainParam.showWindow = 1; 

[net, tr] = train(net, sampl, targ);