

clear all; close all; clc;

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau' ,'b', 'c', 'NR', 'a',...
    'k_tau','k_{c}'};
% "NR" - not relevnt param 

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;

% % file name for uploading:
% results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_All_1.mat'};

%% Load data:
results = load_results(results_fileName);

%% get and filter periods:
% % get groups:
[~,~,~,ids_osc] = get_CPGs(results,'osc','2N',MML);
[~,~,~,ids_n_osc] = get_CPGs(results,'n-osc','2N',MML);

targ = [ids_osc; ids_n_osc];

seq = (vertcat(results(:).seq))';
periods = horzcat(results(:).periods);

% plot histograms:
% plot_osc_nosc_hist([3,3],seq,periods,seqOrder,ids_osc,ids_n_osc);
% plot_param_hist(seq,periods,seqOrder)

% input_names = {'tau','b','c','a'};
% output_names = {'b'};

input_names = {'tau','b','a'};
output_names = {'b'};

% input_names = {'tau','a'};
% output_names = {'b'};

[sampl,~] = ...
    prepare_NN_data(input_names,output_names,...
    seqOrder,seq,periods);

plot_pointCloud_clusters(sampl(1,:),sampl(3,:),sampl(2,:),...
    ids_osc,ids_n_osc,{'tau','a','b'});
%% Neural Network #2:
X=sampl;
T=double(targ);
%Train an autoencoder with a hidden layer of size 10 and a linear transfer function for the decoder. Set the L2 weight regularizer to 0.001, sparsity regularizer to 4 and sparsity proportion to 0.05.

hiddenSize = 3;
autoenc1 = trainAutoencoder(X,hiddenSize,...
    'L2WeightRegularization',0.001,...
    'SparsityRegularization',4,...
    'SparsityProportion',0.05,...
    'DecoderTransferFunction','purelin');
%
%Extract the features in the hidden layer.

features1 = encode(autoenc1,X);
%Train a second autoencoder using the features from the first autoencoder. Do not scale the data.

hiddenSize = 2;
autoenc2 = trainAutoencoder(features1,hiddenSize,...
    'L2WeightRegularization',0.001,...
    'SparsityRegularization',4,...
    'SparsityProportion',0.05,...
    'DecoderTransferFunction','purelin',...
    'ScaleData',false);

features2 = encode(autoenc2,features1);

rand_ids = randsample(1:size(sampl,2),1000);
figure;
gscatter(features2(1,rand_ids),features2(1,rand_ids),...
    {ids_osc(rand_ids),ids_in_range(rand_ids),ids_n_osc(rand_ids)},...
    'bgr','..x');
xlabel('feature1'); ylabel('feature2');
legend('osc','osc in range','n-osc');
clear rand_ids

rand_ids = randsample(1:size(sampl,2),10000);
figure;
gscatter(sampl(1,rand_ids),sampl(2,rand_ids),...
    {ids_osc(rand_ids),ids_in_range(rand_ids),ids_n_osc(rand_ids)},...
    'bgr','..x');
xlabel('tau'); ylabel('b');
legend('osc','osc in range','n-osc');
clear rand_ids

%
softnet = trainSoftmaxLayer(features2,T,'LossFunction','crossentropy');
% %Stack the encoders and the softmax layer to form a deep network.

deepnet = stack(autoenc1,autoenc2,softnet);
%Train the deep network on the wine data.

deepnet = train(deepnet,X,T);
%Estimate the deep network, deepnet.

y = deepnet(X);
plotconfusion(T,deepnet(X));

%% Neural Network #3:
X=sampl;
T=double(targ);
%Train an autoencoder with a hidden layer of size 10 and a linear transfer function for the decoder. Set the L2 weight regularizer to 0.001, sparsity regularizer to 4 and sparsity proportion to 0.05.

hiddenSize = [2];
autoenc1 = trainAutoencoder(X,hiddenSize,...
    'L2WeightRegularization',0.001,...
    'SparsityRegularization',4,...
    'SparsityProportion',0.05,...
    'DecoderTransferFunction','purelin',...
    'ScaleData',false);
%
%Extract the features in the hidden layer.

features1 = encode(autoenc1,X);
%Train a second autoencoder using the features from the first autoencoder. Do not scale the data.

rand_ids = randsample(1:size(sampl,2),1000);
figure;
gscatter(features1(1,rand_ids),features1(1,rand_ids),...
    {ids_osc(rand_ids),ids_n_osc(rand_ids)},...
    'br','.x');
xlabel('feature1'); ylabel('feature2');
legend('osc','n-osc');
clear rand_ids

rand_ids = randsample(1:size(sampl,2),10000);
figure;
gscatter(sampl(1,rand_ids),sampl(2,rand_ids),...
    {ids_osc(rand_ids),ids_n_osc(rand_ids)},...
    'br','.x');
xlabel('tau'); ylabel('b');
legend('osc','n-osc');
clear rand_ids

%
softnet = trainSoftmaxLayer(features1,T,'LossFunction','crossentropy');
% %Stack the encoders and the softmax layer to form a deep network.

deepnet = stack(autoenc1,softnet);
%Train the deep network on the wine data.

deepnet = train(deepnet,X,T);
%Estimate the deep network, deepnet.

y = deepnet(X);
plotconfusion(T,deepnet(X));

%% Neural Network #4: (using 'patternnet')
X=sampl;
T=double(targ);
%Train an autoencoder with a hidden layer of size 10 and a linear transfer function for the decoder. Set the L2 weight regularizer to 0.001, sparsity regularizer to 4 and sparsity proportion to 0.05.

architecture = [10];
net = patternnet(architecture);
[net, tr] = train(net, X, T);

rand_ids = randsample(1:size(sampl,2),10000);
figure;
gscatter(sampl(1,rand_ids),sampl(2,rand_ids),...
    {ids_osc(rand_ids),ids_n_osc(rand_ids)},...
    'br','.x');
xlabel('tau'); ylabel('b');
legend('osc','n-osc');
clear rand_ids

deepnet = net; % for now...

y = deepnet(X);
plotconfusion(T,deepnet(X));

%%
close all;

rand_seq = MML.Gen.RandSeq(1);

[NN_in_test,~] = ...
    prepare_NN_data(input_names,output_names,...
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
[ Seq_old ] = MML.Gen.RandSeq(1000); %generate 'N' rand samples
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
        prepare_NN_data(input_names,output_names,...
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
