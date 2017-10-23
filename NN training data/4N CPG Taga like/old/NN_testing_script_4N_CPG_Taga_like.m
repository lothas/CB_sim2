
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


Mw = 10*ones(1,4);
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

results_fileName = {'MatsRandomRes_TagaLike_CPG_all_1.mat',...
    'MatsRandomRes_TagaLike_CPG_all_2.mat'};
%% Load data:
results_before = load_results(results_fileName);
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
    '    periods1: ',...
    num2str(results(rand_id).periods(1)),...
    '    periods2: ',...
    num2str(results(rand_id).periods(2))]});
subplot(2,1,2)
plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');

clear out signal N rand_id

%% get and filter periods:
% % get oscillating:
[results,periods,seq,~] = get_CPGs(results_before,'osc','2N',MML);
% plot_param_hist(seq,periods,seqOrder)

%% Prepare NN inputs and outputs:
input_names = {'b','tau','a'};
output_names = {'periods'};

[sampl,targ] = ...
    prepare_NN_data(input_names,output_names,...
    seqOrder,seq,periods);

% % % Neural Network #1a: (CPU training)
architecture = [6,20];

net = fitnet(architecture);
% net = feedforwardnet(architecture);
net.trainFcn = 'trainbr';
net.divideParam.trainRatio = 0.5;
net.divideParam.valRatio = 0.35;
net.divideParam.testRatio = 0.15;

net.trainParam.showWindow = 1; 
net.trainParam.showCommandLine = 1;
net.trainParam.epochs = 1000;
% net.performParam.normalization = 'percent'; %It must be 'none', 'standard' or 'percent'

[net, tr] = train(net, sampl, targ);
[r,m,b] = regression(targ,net(sampl));
disp(['r = ',num2str(r),'    m = ',num2str(m),'    b = ',num2str(b)]);
clear r m b

figure;
histogram(targ,100,'Normalization','pdf'); hold on;
histogram(net(sampl),100,'Normalization','pdf');
legend('targets','NN outputs');

%% Classification NN:
NNclasi_in = {'tau','b','c','w_{1}','w_{2}','w_{3}','w_{4}'};
NNclasi_out = {'b'}; % not really relevant

% % get groups:
[~,~,~,ids_osc] = get_CPGs(results_before,'osc','4N',MML);
[~,~,~,ids_n_osc] = get_CPGs(results_before,'n-osc','4N',MML);
targ_clasi = [ids_osc; ids_n_osc];
seq = (vertcat(results_before(:).seq))';
periods = horzcat(results_before(:).periods);

[sampl_clasi,~] = ...
    prepare_NN_data(NNclasi_in,NNclasi_out,...
    seqOrder,seq,periods);

%% Neural Network #4: (using 'patternnet')
X=sampl_clasi;
T=double(targ_clasi);
%Train an autoencoder with a hidden layer of size 10 and a linear transfer function for the decoder. Set the L2 weight regularizer to 0.001, sparsity regularizer to 4 and sparsity proportion to 0.05.

architecture = [7,2];
net = patternnet(architecture);
[net, tr] = train(net, X, T);

rand_ids = randsample(1:size(sampl_clasi,2),10000);
figure;
gscatter(sampl_clasi(1,rand_ids),sampl_clasi(2,rand_ids),...
    {ids_osc(rand_ids),ids_n_osc(rand_ids)},...
    'br','.x');
xlabel('tau'); ylabel('b');
legend('osc','n-osc');
clear rand_ids

deepnet = net; % for now...

y = deepnet(X);
plotconfusion(T,deepnet(X));

%% Test net on random seq:
close all;

rand_seq = MML.Gen.RandSeq(1);

[NN_in_test,~] = ...
    prepare_NN_data(NNclasi_in,NNclasi_out,...
    seqOrder,rand_seq',0.7);

check = deepnet(NN_in_test);
[~,ind] = max(check);
states = {'osc','n-osc'};

[out, ~, signal] = MML.runSim(rand_seq);
figure;
subplot(2,1,1);
plot(signal.T,signal.X);
xlabel('time[sec]');    ylabel('X_i');
title({'X_i over time',...
    ['periods: ',...
    num2str(out.periods'),...
    '     CPG identified: ',states{1,ind}]});
subplot(2,1,2)
plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');
clear signal rand seq out ind states

%% test the NN:
clc;

% % % % CPG parameters:
[ Seq_old ] = MML.Gen.RandSeq(500); %generate 'N' rand samples
% % % run the rand samples to check periods:

% ini the structures to the right size:
clear results_old results_new
results_old(length(Seq_old)).seq = [];
results_old(length(Seq_old)).periods = [];
results_new(length(Seq_old)).seq = [];
results_new(length(Seq_old)).periods = [];

disp('start with the sim:');
parfor i=1:length(Seq_old) % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
% for i=1:length(Seq_old)
    disp(['at sim #',num2str(i)]);
    [out,~, signal] = MML.runSim(Seq_old(i,:));
    
    % Prepare output :
    results_old(i).seq = Seq_old(i,:);
    results_old(i).periods = out.periods;
    
    % don't do anything if CPG IS stable
    if ~any(isnan(out.periods))
        results_new(i).seq = results_old(i).seq;
        results_new(i).periods = out.periods;
        continue;
    end
    
    rand_seq = MML.Gen.RandSeq(1000);
    [NN_in_test,~] = ...
        prepare_NN_data(NNclasi_in,NNclasi_out,...
        seqOrder,rand_seq',0.7);
    check = deepnet(NN_in_test);
    [~,ind] = max(check);
    rand_good_ind = randsample(find(ind==1),1);
    seq_temp = rand_seq(rand_good_ind,:);
    
%     % make NN inputs:
%     rand_seq = Seq_old(i,:);
%     for j=1:100
%         [NN_in_test,~] = ...
%         prepare_NN_data(input_names,output_names,...
%         seqOrder,rand_seq',0.7);
%         
%     % check if CPG is oscillatory:
%         check = deepnet(NN_in_test);
%         [~,ind] = max(check);
%         if ind == 1
%             continue;
%         elseif ind == 3
%             rand_seq = MML.Gen.RandSeq(1);
%         end
%     end
%     seq_temp = rand_seq;
    
    % TODO: add rescale
    
    % run sim again:
    [out, ~, signal] = MML.runSim(seq_temp);
        % Prepare output:
    % Parameters
    results_new(i).seq = seq_temp;

    results_new(i).periods = out.periods;

end 
disp('sim end...');

% get CPG periods:
periods_old = horzcat(results_old(:).periods);
periods_old = periods_old(2,:); % get rid of hip period
% Filter CPG's where not both signals oscillating:
osc_ids_temp = ~isnan(periods_old);

disp(['the num of osc CPGs before the NN: ',num2str(sum(osc_ids_temp))]);

% get CPG periods:
periods_new = horzcat(results_new(:).periods);
periods_new = periods_new(2,:); % get rid of hip period
% Filter CPG's where not both signals oscillating:
osc_ids_temp = ~isnan(periods_new);

disp(['the num of osc CPGs before the NN: ',num2str(sum(osc_ids_temp))]);

seq_old = vertcat(results_old(:).seq);
seq_new = vertcat(results_new(:).seq);

figure;
histogram(seq_old(:,2),100,'Normalization','pdf'); hold on;
histogram(seq_new(:,2),100,'Normalization','pdf');