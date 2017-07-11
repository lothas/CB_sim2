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

periods_AC_mean = mean(periods_AC,1);
periods_LS_mean = mean(periods_LS,1);

% TODO: do something is one of the periods os NaN

%% get varifications:

varification_1_vec = (vertcat(results(:).neuronOsc))';

varification_1 = (varification_1_vec(1,:) | varification_1_vec(2,:)) & ...
    (varification_1_vec(3,:) | varification_1_vec(4,:));

disp(['the number of CPG that pass AC and didnt pass varification #1 is:',...
    ' ',num2str(sum(varification_1 & ~(ids_AC_ankle & ids_AC_hip)))]);
clear varification_1_vec

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

%% the entire consusion matrix:
AC_ids = [ids_AC_hip & ids_AC_ankle;
    ids_AC_hip & ~ids_AC_ankle;
    ~ids_AC_hip & ids_AC_ankle;
    ~ids_AC_hip & ~ids_AC_ankle];

LS_ids = [ids_LS_hip & ids_LS_ankle;
    ids_LS_hip & ~ids_LS_ankle;
    ~ids_LS_hip & ids_LS_ankle;
    ~ids_LS_hip & ~ids_LS_ankle];

figure;
plotconfusion(AC_ids,LS_ids,'total')
xlabel('Auto-Corr')
set(gca,'xticklabel',{'hip osc, ankle osc'...
    'hip osc, ankle n-osc' 'hip n-osc, ankle osc'...
    'hip n-osc, ankle n-osc' ''})
ylabel('LS')
set(gca,'yticklabel',{'hip osc, ankle osc'...
    'hip osc, ankle n-osc' 'hip n-osc, ankle osc'...
    'hip n-osc, ankle n-osc' ''})

%% correlation coefficient between AC and LS
best_ids_ankle = ids_AC_ankle & ids_LS_ankle;
best_ids_hip = ids_AC_hip & ids_LS_hip;

R_ankle = corrcoef(periods_AC_ankle(1,best_ids_ankle),...
    periods_LS_ankle(1,best_ids_ankle))

R_hip = corrcoef(periods_AC_hip(1,best_ids_hip),...
    periods_LS_hip(1,best_ids_hip))
%% repeat with mean periods
ids_AC_mean = ~isnan(periods_AC_mean);
ids_LS_mean = ~isnan(periods_LS_mean);

figure;
plotconfusion(ids_AC_mean,ids_LS_mean,'Mean torque')
xlabel('Auto-Corr')
set(gca,'xticklabel',{'n-osc' 'osc' ''})
ylabel('LS')
set(gca,'yticklabel',{'n-osc' 'osc' ''})

% find samples which AC gives osc and LS gives n-osc:
testIds_1 = find((ids_AC_mean & ~ids_LS_mean));

% find samples which AC gives n-osc and LS gives osc:
testIds_2 = find((~ids_AC_mean & ids_LS_mean));
% plotCPG(results,3274)
%% check conditions from section 3.5 in paper:
y_i_osc = (vertcat(results(:).neuronOsc))';

y_12_osc = y_i_osc(1,:) | y_i_osc(2,:);
y_34_osc = y_i_osc(3,:) | y_i_osc(4,:);

y_12_and_34_osc = y_12_osc & y_34_osc;

periods_AC = horzcat(results(:).periods);
osc_net = ~isnan(periods_AC(1,:)) & ~isnan(periods_AC(2,:));

%% plot example:
close all
plotCPG(results,'random')

%% check Matsuoka conditions:
clc
nSims = length(results);

cond0 = false(1,nSims);
cond1 = false(4,nSims);
cond2 = false(1,nSims);

for i=1:nSims
    
    cond0(1,i) = checkCond_0(results(i).seq,seqOrder);
    
    cond1_temp = checkCond_1(results(i).seq,seqOrder);
    cond1(:,i) = cond1_temp';
    
    cond2(1,i) = checkCond_2(results(i).seq,seqOrder,cond1_temp);
    
    if mod(i,10000)==0
        disp(['at i=',num2str(i),...
            '  cond0=',num2str(cond0(1,i)),...
            '  cond1=[',num2str((cond1(:,i))'),']',...
            '  cond2=',num2str(cond2(1,i))]);
    end
end

disp(['CPGs which fulfil cond #0 : ',num2str(sum(cond0))]);

disp(['CPGs which fulfil cond #2 : ',num2str(sum(cond2))]);

% taking only samples with both non-NaN periods:
ids_AC_mean = ~isnan(mean(periods_AC,1));

% plot confusion matrix for condition2:
figure;
plotconfusion(cond2,ids_AC_mean,'condition 2')
xlabel('cond2')
set(gca,'xticklabel',{'false' 'true' ''})
ylabel('CPG osc')
set(gca,'yticklabel',{'false' 'true' ''})
