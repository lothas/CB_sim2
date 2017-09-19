
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

% file name for uploading:
results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_in_des_range_1.mat'};

% results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_only_osc_1.mat'};

% results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_smaller_per_range_1mat'};

% results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W.mat'};

% results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_uniform_b_1.mat',...
%     'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_uniform_b_2.mat'};
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
% % get oscillating:
[results,periods,seq,~] = get_CPGs(results_before,'osc',MML);
% 
% % % also filter out abnormaly big periods:
% good_ids = periods < 5;
% results = results(good_ids);
% periods = periods(:,good_ids);
% seq = seq(:,good_ids);
% clear good_ids

% % % get oscillating in period range:
% [results,periods,seq,~] = get_CPGs(results_before,'osc_in_per_range',MML);

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

% % % normalizing the inputs and outputs:
% [sampl_n,sampl_s] = mapstd(sampl);
% [targ_n,targ_s] = mapstd(targ);

%% Neural Network #1:
architecture = [20];

net = fitnet(architecture);
% net = feedforwardnet(architecture);
net.trainFcn = 'trainbr';
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

% [net, tr] = train(net, sampl_n, targ_n);
% figure;
% histogram(targ_n,100,'Normalization','pdf'); hold on;
% histogram(net(sampl_n),100,'Normalization','pdf');
% legend('targets','NN outputs');
%% Neural Network #2:
X=sampl ;
T=targ;
%Train an autoencoder with a hidden layer of size 10 and a linear transfer function for the decoder. Set the L2 weight regularizer to 0.001, sparsity regularizer to 4 and sparsity proportion to 0.05.

hiddenSize = 2;
autoenc1 = trainAutoencoder(X,hiddenSize,...
    'L2WeightRegularization',0.001,...
    'SparsityRegularization',4,...
    'SparsityProportion',0.05,...
    'DecoderTransferFunction','purelin');
%
%Extract the features in the hidden layer.

features1 = encode(autoenc1,X);
%Train a second autoencoder using the features from the first autoencoder. Do not scale the data.

hiddenSize = 1;
autoenc2 = trainAutoencoder(features1,hiddenSize,...
    'L2WeightRegularization',0.001,...
    'SparsityRegularization',4,...
    'SparsityProportion',0.05,...
    'DecoderTransferFunction','purelin',...
    'ScaleData',false);

features2 = encode(autoenc2,features1);

%
% softnet = trainSoftmaxLayer(features2,T,'LossFunction','crossentropy');
% %Stack the encoders and the softmax layer to form a deep network.
net1 = fitnet(10);
net1.trainFcn = 'trainbr';
net1 = train(net1,features2,T);

deepnet = stack(autoenc1,autoenc2,net1);
%Train the deep network on the wine data.

deepnet = train(deepnet,X,T);
%Estimate the deep network, deepnet.

y = deepnet(X);
net = deepnet;
% clear X T y deepnet softnet features2 autoenc2 hiddenSize features1 autoenc1
%%
clc;% close all
% 
% input_names = {'tau','a','periods'};
% output_names = {'b'};
% plot_param_hist(seq,periods,seqOrder)

N = 100;
M = 5;

tau = 0.03;
a = linspace(1,3.5,M);
p = linspace(0.1,15,N);%0.7;

NNout = zeros(M,N);

for j=1:M
    for i=1:N
        NNin = [tau;a(1,j);p(1,i)];
        NNout(j,i) = net(NNin);
%         NNin_n = mapstd('apply',NNin,sampl_s);
%         NNout_n = sim(net,NNin_n);
%         NNout(j,i) = mapstd('reverse',NNout_n,targ_s);        
    end
end

figure; hold on;
for j=1:M
    plot(p,NNout(j,:),'Marker','o');
    LABELS{1,j} = sprintf('a = %f',a(1,j));
end
grid minor;
legend(LABELS);
title(['tau = ',num2str(tau)]);
xlabel('periods');
ylabel('NN outputs');

clear N tau a p NNout LABELS rand_id N
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

% % % % % find nearest neighbours:
[IDX,D] = knnsearch(NN_out',targ');

figure;
scatter(NN_out(1:1000),targ(1,IDX(1:1000)));
xlabel('NN output "b"');
ylabel('nearest neighbor from the training data')
% axis([0 2.5 0 2.5])
% % % % % % % % % % % % % % % 

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