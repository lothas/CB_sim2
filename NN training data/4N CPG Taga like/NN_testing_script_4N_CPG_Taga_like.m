
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

% % % % % % % % % Narrow b, Narrow W
% Mw = 5*ones(1,4);
%     % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
% Keys = {'\tau_r', 'beta', 'amp',   '4neuron_taga_like', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%               1 ,      1,    4 ,           4,        1 ,         4 ,            0 };
% Range = {  0.02 ,    0.2,  mamp,          mw,      -10 ,  -0.1*Mamp; % Min
%            0.25 ,    2.5,  Mamp,          Mw,       10 ,   0.1*Mamp}; % Max

% % % % % % % % % Narrow b, Large W
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

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau','b','c_1','c_2','c_3','c_4',... 
'w_{1}','w_{2}','w_{3}','w_{4}'};

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;

% % change tau_a/tau_r to 12 (instead of 5)
MML.Sim.Con.tau_ratio = 12;

% file name for uploading:L
% results_fileName = {'MatsRandomRes_4Neurons_TagaLike_Narrow_b_Narrow_W_1.mat'};
results_fileName = {'MatsRandomRes_4Neurons_TagaLike_Narrow_b_Large_W_1.mat'};
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
    num2str(results(rand_id).periods(1)),...
    ' ',num2str(results(rand_id).periods(2))]});
subplot(2,1,2)
plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');

clear out signal N rand_id

%% get and filter periods:

% get CPG periods:
periods = horzcat(results(:).periods);

% Filter CPG's where not both signals oscillating:
osc_ids_temp = ~isnan(periods);
osc_ids_temp = osc_ids_temp(1,:) & osc_ids_temp(2,:);
disp(['Number of non-osc CPGs: ',num2str(sum(~osc_ids_temp))]);

% Filter CPG's where the is a big difference between hip and ankle:
periods_ratios = (periods(1,:)./periods(2,:));
diff_ids = (periods_ratios >  0.85) & (periods_ratios <  1.15); 
disp(['Number of CPGs with not matching periods: (from the osc ones)',...
    num2str(sum(osc_ids_temp & ~diff_ids))]);
periods = mean(periods,1);

% % plot the distribution of the missdefined CPG periods:
if false
    figure;
    h=histogram(periods_ratios,100); grid minor;
    h.BinLimits = [0,2.5];
    h.BinWidth = 0.1;
    h.Normalization = 'pdf';
    xlabel('P_{hip} / P_{ankle} : the ratio between the two CPG outputs');
    ylabel('probability density');
    title('histogram of the ratio between the two CPG outputs');
    set(gca,'FontSize',10);
    savefig('figure_TBD_Histogram_of_ratio_between_periods_hipAnkle')
end

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

num_of_osc_ids_exluded_param_range = sum(osc_ids_temp & diff_ids);
osc_ids = osc_ids_temp & diff_ids & ids_in_genome_range;

osc_inRange_ids = osc_ids &...
    ( (periods(1,:) > MML.perLimOut(1,1)) &...
    (periods(1,:) < MML.perLimOut(1,2)) );

num_of_osc_ids = sum(osc_ids);
num_of_inRange_ids = sum(osc_inRange_ids);

%% keep the good seq and periods
seq = seq(:,osc_ids);
periods = periods(:,osc_ids);

%% Prepare NN inputs and outputs:
input_names = {'periods','tau',...
    'w_{1}','w_{2}','w_{3}','w_{4}'};

output_names = {'b'};

[sampl,targ] = ...
    prepare_NN_data(input_names,output_names,...
    seqOrder,seq,periods);

%% Neural Network:
architecture = [10];

net = fitnet(architecture);
net.trainFcn = 'trainbr';
net.trainParam.showWindow = 1; 

[net, tr] = train(net, sampl, targ);