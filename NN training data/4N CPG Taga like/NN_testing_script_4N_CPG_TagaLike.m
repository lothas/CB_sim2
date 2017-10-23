
%% NN testing Script:
% IMPORTANT: dont forget to load the right genome file and to uptade
% 'MatsuokaML.m' to the right rettings
%
% Seriously: Don't forget to check in each fileheader the right gene
% boundaries.

clear all; close all; clc;

% % % the order of the parametrs in CPG Sequence:
% Symm tonic inputs:
seqOrder = {'tau' ,'b', 'c', 'w1', 'w2', 'w3', 'w4',...
    'k_tau','k_{c}'};
% % general tonic inputs:
% seqOrder = {'tau' ,'b', 'c1', 'c2', 'c3', 'c4',...
%     'w1','w2','w3','w4',...
%     'k_tau','k_{c1}','k_{c2}','k_{c3}','k_{c4}'};

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;

results_fileName = {'MatsRandomRes_TagaLike_TrainingSet.mat'};
%% Load data:
results = load_results(results_fileName);

%% plot random example:
clc; close all

N = length(results);

% % show random CPG:
rand_id = randsample(1:N,1);
[out, ~, signal] = MML.runSim(results(rand_id).seq);

% figure;
% subplot(2,1,1);
% plot(signal.T,signal.X);
% xlabel('time[sec]');    ylabel('X_i');
% title({'X_i over time',...
%     ['id #',num2str(rand_id),...
%     '    periods: ',...
%     num2str(results(rand_id).periods(2))]});
% subplot(2,1,2)
% plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');
% 
% clear out signal rand_id

[~,~,~,ids_osc] = get_CPGs(results,'osc','4N',MML);

figure;
for i=1:25
%     rand_id = randsample(1:N,1);
    rand_id = randsample(find(ids_osc),1);
    [~, ~, signal] = MML.runSim(results(rand_id).seq);
    subplot(5,5,i);
    plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');
    xlabel('time[sec]');    ylabel('CPG output');
    title( ['id #',num2str(rand_id),...
        '    periods: ',...
        num2str(results(rand_id).periods(1)),'    ',...
        num2str(results(rand_id).periods(2))]);
end
clear i signal rand_id N 
clear ids_osc
%% Regression NN:
% % get oscillating:
[results,periods,seq,~] = get_CPGs(results,'osc','2N',MML);

% 
% figure;
% boxplot(seq','orientation','horizontal','labels',seqOrder)
% 
% plot_param_hist(seq,periods,seqOrder)

% % Prepare NN inputs and outputs:
input_names = {'b','tau','w1','w2'};
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

% [net, tr] = train(net, sampl_n, targ_n);
% figure;
% histogram(targ_n,100,'Normalization','pdf'); hold on;
% histogram(net(sampl_n),100,'Normalization','pdf');
% legend('targets','NN outputs');

%% Classification NN:
% % get groups:
[~,~,~,ids_osc] = get_CPGs(results,'osc','4N',MML);
[~,~,~,ids_n_osc] = get_CPGs(results,'n-osc','4N',MML);

targ = [ids_osc; ids_n_osc];

seq = (vertcat(results(:).seq))';
periods = horzcat(results(:).periods);

% plot histograms:
% plot_osc_nosc_hist([3,3],seq,periods,seqOrder,ids_osc,ids_n_osc);
% plot_param_hist(seq,periods,seqOrder)

input_names = {'tau','b','w1','w2','w3','w4'};
output_names = {'b'};


[sampl,~] = ...
    prepare_NN_data(input_names,output_names,...
    seqOrder,seq,periods);

% % define the NN:
X=sampl;
T=double(targ);
architecture = [10];
net = patternnet(architecture);
[net, tr] = train(net, X, T);

deepnet = net; % for now...

y = deepnet(X);
plotconfusion(T,deepnet(X));

%% Check classification net prediction:
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

%% test classification NN [erformance:
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
    
    rand_seq = MML.Gen.RandSeq(10000);
    [NN_in_test,~] = ...
        prepare_NN_data(input_names,output_names,...
        seqOrder,rand_seq',0.7);
    check = deepnet(NN_in_test);
    [~,ind] = max(check);
    rand_good_ind = randsample(find(ind==1),1);
    seq_temp = rand_seq(rand_good_ind,:);
    
    % run sim again:
    [out, ~, signal] = MML.runSim(seq_temp);
        % Prepare output:
    % Parameters
    results_new(i).seq = seq_temp;

    results_new(i).periods = out.periods;

end 
disp('sim end...');

% get CPG periods:
[~,~,~,ids_osc_old] = get_CPGs(results_old,'osc','4N',MML);
disp(['the num of osc CPGs before the NN: ',num2str(sum(ids_osc_old))]);

% get CPG periods:
[~,~,~,ids_osc_new] = get_CPGs(results_new,'osc','4N',MML);
disp(['the num of osc CPGs after the NN: ',num2str(sum(ids_osc_new))]);

%% plot random example (after NN):
clc; close all

N = length(results_new);

figure;
for i=1:25
    rand_id = randsample(1:N,1);
    [~, ~, signal] = MML.runSim(results_new(rand_id).seq);
    subplot(5,5,i);
    plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');
    xlabel('time[sec]');    ylabel('CPG output');
    title( ['id #',num2str(rand_id),...
        '    periods: ',...
        num2str(results_new(rand_id).periods(1)),'    ',...
        num2str(results_new(rand_id).periods(2))]);
end
clear i signal rand_id N 
clear ids_osc