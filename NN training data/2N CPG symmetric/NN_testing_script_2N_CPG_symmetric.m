
%% NN testing Script:
% IMPORTANT: dont forget to load the right genome file and to uptade
% 'MatsuokaML.m' to the right rettings
% 

clear all; close all; clc;

%% Create genome (only if necessary)
genome_file = 'MatsuokaGenome_2Neuron_Symm.mat';
nAnkle = 1;%1; % Number of ankle torques
nHip = 0;   % Number of hip torques
maxAnkle = 10;   % Max ankle torque
maxHip = 10;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;

% %     %rrr 2neuron symmetric specific range%%
% % large b Large W
% Mw = 10;
% mw = 0;
% Keys = {'\tau_r', 'beta',     'amp_2n',    '2neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%               1 ,      1,            2,                         1,        1 ,          2,            0 };
% Range = {  0.02 ,    0.2,        [0,0],                         0,   -0.001 ,[-0.2,-0.2]; % Min
%            0.25  ,    10,      [10,10],                        10,    0.001 , [0.2,0.2]}; % Max

%        % Narrow b Large W
% Mw = 10;
% mw = 0;
% Keys = {'\tau_r', 'beta',     'amp_2n',    '2neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
%               1 ,      1,            2,                         1,        1 ,          2,            0 };
% Range = {  0.02 ,    0.2,        [0,0],                         0,   -0.001 ,[-0.2,-0.2]; % Min
%            0.25  ,   2.5,      [10,10],                        10,    0.001 , [0.2,0.2]}; % Max

       % Narrow b narrow W
Mw = 5;
mw = 0;
Keys = {'\tau_r', 'beta',     'amp_2n',    '2neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
              1 ,      1,            2,                         1,        1 ,          2,            0 };
Range = {  0.02 ,    0.2,        [0,0],                         0,   -0.001 ,[-0.2,-0.2]; % Min
           0.25  ,   2.5,      [10,10],                         5,    0.001 , [0.2,0.2]}; % Max

MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

clear all
%%

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau' ,'b', 'c', 'NR', 'a',...
    'k_tau','k_{c1}','k_{c2}'};
% "NR" - not relevnt param 

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 2;

% fix a problem with the not relevant parameter in the seq:
MML.Gen.Range(1,4) = -1;

% file name for uploading:L
% results_fileName = {'MatsRandomRes_2Neurons_symm_Large_b_Large_W.mat'};
% results_fileName = {'MatsRandomRes_2Neurons_symm_4Paper.mat'};
results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W.mat'};
%% Load data:
load(results_fileName{1,1},'results','header');
disp('data file information:');
disp(header);
% if I have more than 1 data file:
for i=2:numel(results_fileName)
    data = load(results_fileName{1,i},'results');
    results = [results, data.results]; %#ok<AGROW>
end

results_before = results;
clear results
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
% % % get oscillating:
% [results,periods,seq] = get_CPGs(results_before,'osc',MML);
% 
% % % also filter out abnormaly big periods:
% good_ids = periods < 5;
% results = results(good_ids);
% periods = periods(:,good_ids);
% seq = seq(:,good_ids);
% clear good_ids

% get oscillating in period range:
[results,periods,seq] = get_CPGs(results_before,'osc_in_per_range',MML);

figure;
boxplot(seq','orientation','horizontal','labels',seqOrder)

plot_param_hist(seq,periods,seqOrder)

%% remove outliers:
% https://www.mathworks.com/matlabcentral/answers/121247-how-can-i-detect-and-remove-outliers-from-a-large-dataset

ids = true(1,size(seq,2));

for i=1:(MML.Gen.Length+1)
    
    if i > MML.Gen.Length
        vector = periods;
    else
        vector = seq(i,:);
    end

    percntiles = prctile(vector,[5 95]); %5th and 95th percentile
    % the distance between the %5th and 95th percentiles is four stdevs
    
    outlierIndexes = vector < percntiles(1) | vector > percntiles(2);
    
    ids = ~outlierIndexes & ids;
%     % Extract outlier values:
%     outliers = vector(outlierIndexes);
%     % Extract non-outlier values:
%     nonOutliers = vector(~outlierIndexes);
end

periods = periods(:,ids);
seq = seq(:,ids);

figure;
subplot(2,1,1);
boxplot(seq','orientation','horizontal','labels',seqOrder);
subplot(2,1,2);
boxplot(periods,'orientation','horizontal','labels',{'periods'});

clear ids
%% Prepare NN inputs and outputs:
% input_names = {'b','tau','a'};
% output_names = {'periods'};

input_names = {'tau','a','periods'};
output_names = {'b'};

[sampl,targ] = ...
    prepare_NN_data(input_names,output_names,...
    seqOrder,seq,periods);

%% Neural Network:
architecture = [5];

net = fitnet(architecture);
net.trainFcn = 'trainbr';
net.trainParam.showWindow = 1; 

[net, tr] = train(net, sampl, targ);

figure;
histogram(targ,100,'Normalization','pdf'); hold on;
histogram(net(sampl),100,'Normalization','pdf');
legend('targets','NN outputs');

%%
clc; close all

N = 100;
M = 5;

tau = 0.03;
a = linspace(0.1,5,M);
p = linspace(0.68,0.78,N);%0.7;

NNout = zeros(M,N);

for j=1:M
    for i=1:N
        NNout(j,i) = net([tau;a(1,j);p(1,i)]);
    end
end
N = length(results);
rand_id = randsample(1:N,1);

% % seqOrder = {'tau' ,'b', 'c', 'NR', 'a',...
% %     'k_tau','k_{c1}','k_{c2}'};
% [~, ~, signal] = MML.runSim([0.021,2,7,0.001,2,0,0,0]);
% figure;
% subplot(2,1,1);
% plot(signal.T,signal.X);
% xlabel('time[sec]');    ylabel('X_i');
% subplot(2,1,2)
% plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');
% clear signal

figure; hold on;
for j=1:M
    plot(p,NNout(j,:),'Marker','o');
    LABELS{1,j} = sprintf('a = %f',a(1,j));
end
grid minor;
legend(LABELS);
xlabel('periods');
ylabel('NN outputs');

clear N tau a p NNout LABELS
%% % test NN with training data (only with des period)
close all

desPeriod = MML.perLim(1) + ...
                 rand(1,size(targ,2))*(MML.perLim(2)-MML.perLim(1));
[NN_in,~] = ...
    prepare_NN_data(input_names,output_names,...
    seqOrder,seq,desPeriod);

% % % %  test on rande seq #1:
% rand_seq = MML.Gen.RandSeq(length(targ));
% [NN_in,~] = ...
%     prepare_NN_data(input_names,output_names,...
%     seqOrder,rand_seq',desPeriod);
% figure;
% boxplot(rand_seq,'orientation','horizontal','labels',seqOrder)
% clear rand_seq

NN_out = net(NN_in);
plotregression(targ,NN_out)

% NN_out = net(sampl(:,tr.testInd));
% plotregression(targ(:,tr.testInd),NN_out)

figure;
histogram(NN_out,100,'Normalization','pdf');
title('histogram of NN_{output}');
xlabel('NN_{output}');

figure;
histogram(periods,100,'Normalization','pdf');
title('histogram of periods');
xlabel('periods');

figure;
subplot(2,1,1);
boxplot(sampl','orientation','horizontal','labels',input_names)
subplot(2,1,2);
boxplot(targ','orientation','horizontal','labels',output_names)

clear NN_in NN_out desPeriod