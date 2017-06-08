
close all; clear all; clc

% filename = 'MatsRandomRes_2Neurons_symm_test_for_FFT.mat';
filename = 'MatsRandomRes_2Neurons_symm_test_for_FFT_only_conv.mat';
load(filename,'results');

% taking only CPGs with oscilations: 
periods = vertcat(results(:).periods);
ids_Jon = (~(isnan(periods)))';
ids_LSQ = false(1,length(periods));
for i=1:length(periods)
    temp = results(i).periods_LSQ;
    ids_LSQ(1,i) = ~(isnan(temp(1,1)));
end
ids_good_period = ids_Jon & ids_LSQ;

% taking only CPGs with reasonable coeffients:
ids_good_bias = true(1,length(periods));
ids_good_Sine = true(1,length(periods));
ids_good_cos = true(1,length(periods));
for i=1:length(periods)
    % check bias coef:
    temp = results(i).bias_coef;
    if abs(temp) > 10
        ids_good_bias(1,i) = false;
    end
    
    % check Sine coef:
    temp = results(i).sine_coef;
    if (abs(temp(1,1)) > 10) || (abs(temp(1,2)) > 10) || (abs(temp(1,2)) > 10)
        ids_good_Sine(1,i) = false;
    end
    
    % check cosSine coef:
    temp = results(i).cos_coef;
    if (abs(temp(1,1)) > 10) || (abs(temp(1,2)) > 10) || (abs(temp(1,2)) > 10)
        ids_good_cos(1,i) = false;
    end
end
ids_good_coef = ids_good_bias & ids_good_Sine & ids_good_cos;

ids = ids_good_period & ids_good_coef;

periods_LSQ = vertcat(results(ids).periods_LSQ);

tau = vertcat(results(ids).tau);
b = vertcat(results(ids).b);
c = vertcat(results(ids).c);
a = vertcat(results(ids).a);


freq_1st = 1./periods_LSQ(:,1); % because the other harmonics are just multiplication
bias_coef = vertcat(results(ids).bias_coef);
sine_coef = vertcat(results(ids).sine_coef);
cos_coef = vertcat(results(ids).cos_coef);

figure;
subplot(2,2,1)
hist(bias_coef,100); grid minor;
title('bias coef hist');
subplot(2,2,2)
hist(sine_coef(:,1),100); grid minor;
title('Sine coef hist');
subplot(2,2,3)
hist(cos_coef(:,1),100); grid minor;
title('Cosine coef hist');


%% Train NN with 1 layer
% inputs = ([tau,b,c,a])';
% targets = ([freq_1st,bias_coef,sine_coef(:,1),cos_coef(:,1)])';

inputs = ([freq_1st,bias_coef,sine_coef(:,1),cos_coef(:,1),b])';
targets = ([tau,c,a])';

HiddenN = 10;
net = feedforwardnet(HiddenN);
% net = fitnet(HiddenN);
[net, tr] = train(net, inputs, targets);


%% try to reconstruct a signal:
freq_des = 0.7;
bias_des = 0;
sine_coef1_des = 1.3;
cos_coef1_des = 1.8;
b_have = 4;
wanted_param = [freq_des;...
    bias_des;...
    sine_coef1_des;...
    cos_coef1_des;...
    b_have];

necessary_out = net(wanted_param);
tau_need = necessary_out(1,1);
c_need = necessary_out(2,1);
a_have = necessary_out(3,1);
% define Matsuoka Sim the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 50; % 15
MML.nNeurons = 2;

seq = [tau_need,b_have,c_need,0,a_have,0,0,0];
[out, ~, signal] = MML.runSim(seq);

x = signal.signal(1,:);
t = signal.T;

figure;
plot(t,x); grid minor;

[~,~,~,~,fitObject] = ...
    MML.processResults_LSQ(x,t,1)
%% Script to unite different samples files

% close all; clear all; clc
% results1 = [];
% for i=1:5
%     filename = ['MatsRandomRes_2Neurons_symm_test_for_FFT_',...
%         num2str(i),'.mat'];
%     load(filename,'results');
%     results1 = [results1,results];
% end
% clear results
% 
% periods = vertcat(results1(:).periods);
% ids_Jon = (~(isnan(periods)))';
% 
% ids_LSQ = false(1,length(periods));
% for i=1:length(periods)
%     temp = results1(i).periods_LSQ;
%     ids_LSQ(1,i) = ~(isnan(temp(1,1)));
% end
% ids = ids_Jon & ids_LSQ;
% 
% results = results1(ids);
% 
% save('MatsRandomRes_2Neurons_symm_test_for_FFT_only_conv.mat','results');