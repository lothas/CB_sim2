
close all; clear all; clc

% filename = 'MatsRandomRes_2Neurons_symm_test_for_FFT.mat';
filename = 'MatsRandomRes_2Neurons_symm_test_for_FFT_only_conv.mat';
load(filename,'results');

periods = vertcat(results(:).periods);
ids_Jon = (~(isnan(periods)))';

ids_LSQ = false(1,length(periods));
for i=1:length(periods)
    temp = results(i).periods_LSQ;
    ids_LSQ(1,i) = ~(isnan(temp(1,1)));
end
ids = ids_Jon & ids_LSQ;

periods_LSQ = vertcat(results(ids).periods_LSQ);

tau = vertcat(results(ids).tau);
b = vertcat(results(ids).b);
c = vertcat(results(ids).c);
a = vertcat(results(ids).a);


freq_1st = 1./periods_LSQ(:,1); % because the other harmonics are just multiplication
bias_coef = vertcat(results(ids).bias_coef);
sine_coef = vertcat(results(ids).sine_coef);
cos_coef = vertcat(results(ids).cos_coef);



%% Train NN with 1 layer
% inputs = ([tau,b,c,a])';
% targets = ([freq_1st,bias_coef,sine_coef(:,1),cos_coef(:,1)])';

inputs = ([freq_1st,bias_coef,sine_coef(:,1),cos_coef(:,1),b,a])';
targets = ([tau,c])';

HiddenN = 3;
net = feedforwardnet(HiddenN);
% net = fitnet(HiddenN);
[net, tr] = train(net, inputs, targets);

% try to reconstruct a signal:
freq_des = 0.7;
bias_des = 0;
sine_coef1_des = 1.3;
cos_coef1_des = 0.64;
b_des = 4;
a_des = 1.8;
wanted_param = [freq_des;...
    bias_des;...
    sine_coef1_des;...
    cos_coef1_des;...
    b_des;...
    a_des];

necessary_out = net(wanted_param);
tau_need = necessary_out(1,1);
c_need = necessary_out(2,1);

% define Matsuoka Sim the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 50; % 15
MML.nNeurons = 2;

seq = [tau_need,b_des,c_need,0,a_des,0,0,0];
[out, ~, signal] = MML.runSim(seq);

x = signal.signal(1,:);
t = signal.T;

figure;
plot(t,x); grid minor;

[~,~,~,~,fitObject] = ...
    MML.processResults_LSQ(x,t,1)
%% Script to unite different samples files

close all; clear all; clc
results1 = [];
for i=1:5
    filename = ['MatsRandomRes_2Neurons_symm_test_for_FFT_',...
        num2str(i),'.mat'];
    load(filename,'results');
    results1 = [results1,results];
end
clear results

periods = vertcat(results1(:).periods);
ids_Jon = (~(isnan(periods)))';

ids_LSQ = false(1,length(periods));
for i=1:length(periods)
    temp = results1(i).periods_LSQ;
    ids_LSQ(1,i) = ~(isnan(temp(1,1)));
end
ids = ids_Jon & ids_LSQ;

results = results1(ids);

save('MatsRandomRes_2Neurons_symm_test_for_FFT_only_conv.mat','results');