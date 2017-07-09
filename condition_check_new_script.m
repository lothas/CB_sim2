clear all; close all; clc

%% Load Data:
results1 = load('MatsRandomRes_4Neurons_with_LSQ_1.mat','results');
results2 = load('MatsRandomRes_4Neurons_with_LSQ_2.mat','results');

results = horzcat(results1.results,results2.results);
clear results1 results2

%% get periods:
periods_AC = horzcat(results(:).periods);
periods_LS = (vertcat(results(:).periods_LSQ))';

periods_AC_hip = periods_AC(1,:);
periods_AC_ankle = periods_AC(2,:);

periods_LS_hip = periods_LS(1,:);
periods_LS_ankle = periods_LS(2,:);

period_AC_mean = mean(periods_AC,1);
period_LS_mean = mean(periods_LS,1);
% TODO: do something is one of the periods os NaN
%% phase 1- confusion matrix, osc/non-osc:
% for ankle torque:
ids_AC_ankle = ~isnan(periods_AC_ankle);
ids_LS_ankle = ~isnan(periods_LS_ankle);

figure;
plotconfusion(ids_AC_ankle,ids_LS_ankle,'Ankle torque')
xlabel('Auto-Corr')
set(gca,'xticklabel',{'n-osc' 'osc' ''})
ylabel('LS')
set(gca,'yticklabel',{'n-osc' 'osc' ''})

% for ankle torque:
ids_AC_hip = ~isnan(periods_AC_hip);
ids_LS_hip = ~isnan(periods_LS_hip);

figure;
plotconfusion(ids_AC_hip,ids_LS_hip,'Hip torque')
xlabel('Auto-Corr')
set(gca,'xticklabel',{'n-osc' 'osc' ''})
ylabel('LS')
set(gca,'yticklabel',{'n-osc' 'osc' ''})

%% correlation coefficient between AC and LS
% ? ? ?  ? ?

%% repeat with mean periods
ids_AC_mean = ~isnan(periods_AC_mean);
ids_LS_mean = ~isnan(periods_LS_mean);

figure;
plotconfusion(ids_AC_mean,ids_LS_mean,'Mean torque')
xlabel('Auto-Corr')
set(gca,'xticklabel',{'n-osc' 'osc' ''})
ylabel('LS')
set(gca,'yticklabel',{'n-osc' 'osc' ''})