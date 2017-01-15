clear all; close all; clc


%% Load Data
load('MatsRandomRes_2Neurons_10_01_2017_A.mat','results')
results1 = results; periods1 = horzcat(results(:).periods);
clear results
load('MatsRandomRes_2Neurons_10_01_2017_B.mat','results')
results2 = results; periods2 = horzcat(results(:).periods);
clear results
load('MatsRandomRes_2Neurons_10_01_2017_C.mat','results')
results3 = results; periods3 = horzcat(results(:).periods);
clear results 
load('MatsRandomRes_2Neurons_10_01_2017_D.mat','results')
results4 = results; periods4 = horzcat(results(:).periods);
clear results 
load('MatsRandomRes_2Neurons_10_01_2017_E.mat','results')
results5 = results; periods5 = horzcat(results(:).periods);
clear results 
load('MatsRandomRes_2Neurons_10_01_2017_F.mat','results')
results6 = results; periods6 = horzcat(results(:).periods);
clear results 
load('MatsRandomRes_2Neurons_10_01_2017_G.mat','results')
results7 = results; periods7 = horzcat(results(:).periods);
clear results 
load('MatsRandomRes_2Neurons_10_01_2017_H.mat','results')
results8 = results; periods8 = horzcat(results(:).periods);
clear results 
load('MatsRandomRes_2Neurons_10_01_2017_I.mat','results')
results9 = results; periods9 = horzcat(results(:).periods);
clear results 
load('MatsRandomRes_2Neurons_10_01_2017_J.mat','results')
results10 = results; periods10 = horzcat(results(:).periods);
clear results 
load('MatsRandomRes_2Neurons_10_01_2017_K.mat','results')
results11 = results; periods11 = horzcat(results(:).periods);
clear results 
load('MatsRandomRes_2Neurons_10_01_2017_L.mat','results')
results12 = results; periods12 = horzcat(results(:).periods);
clear results 
load('MatsRandomRes_2Neurons_10_01_2017_M.mat','results')
results13 = results; periods13 = horzcat(results(:).periods);
clear results 

% concetrate the data:
results = horzcat(results1,results2,results3,results4,results5,results6,...
    results7,results8,results9,results10,results11,results12,results13);
periods = horzcat(periods1,periods2,periods3,periods4,periods5,periods6,...
    periods7,periods8,periods9,periods10,periods11,periods12,periods13);

clear results1 results2 results3 results4 results5 results6 results7...
    results8 results9 results10 results11 results12 results13
clear periods1 periods2 periods3 periods4 periods5 periods6 periods7...
    periods8 periods9 periods10 periods11 periods12 periods13

%% identify which neurons have periods
ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);

%% prepare NN inputs and outputs:
parametersCells = {'tau','b','a','s'};
% parametersCells = {'tau','b'};
targetCells = {'freq'};
[ sampl,targ ] = prepareData_2Neurons( results,periods,ids,parametersCells,targetCells);
%
figure;
histogram(max(horzcat(results(ids).perError2)',[],2),1000);
title('Hist of period error'); xlabel('period error');

clear ids_period ids_error nSims
% clear results periods

%% save to xlse file for eureqa:
filename = 'dataMatrixXLSX_2neurons_updated.xlsx';
dataMatrix = [sampl;targ];
dataMatrix = dataMatrix(:,1:10000);
xlswrite(filename,dataMatrix');
%% normalized inputs
normParams = zeros(size(sampl, 1), 2);
notNormSampl = sampl;
for i = 1:size(sampl, 1)
    feat = sampl(i, :);
    normParams(i, :) = [mean(feat), std(feat)];
    sampl(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
    
end

%% Performance over Hidden neuron Number
NumOfRepeats = 10;
HiddenN = [2,3,4,5,6,7,8,9,10,11,12];
[ NN_Mean_over_HN_num,NN_stdev_over_HN_num ] = NN_Perf_over_HNnum( NumOfRepeats,HiddenN,sampl,targ,1 );

%% Performance over samples quantity
NumOfRepeats = 5;
HiddenN = 5;
dataPointsNum = [1000,5000,7000,10000,15000,20000,30000,50000,70000,...
    100000,120000,150000,164000];
[ NN_Mean_over_sampl_num,NN_stdev_over_sampl_num ] = NN_Perf_over_sampl_num( NumOfRepeats,HiddenN,dataPointsNum,sampl,targ,1 );


%% visualiesd the results in the W matrix
NN_weights_matrix_plot( net,parametersCells )

%% Train NN with 1 layer
HiddenN = 5;
net = feedforwardnet(HiddenN);
% net = fitnet(HiddenN);
[net, tr] = train(net, sampl, targ);

%% Train NN with 2 layers
HiddenN = 10;
net = feedforwardnet([HiddenN,HiddenN]);
[net, tr] = train(net, sampl, targ);
%% Train NN with 3 layers
HiddenN = 3;
net = feedforwardnet([HiddenN,HiddenN,HiddenN]);
[net, tr] = train(net, sampl, targ);

%% Calculating the R^2:

inputs = sampl;
outputs = targ;
netOut = net(inputs);
err = outputs-netOut;
errVar = var(err,0,2);
inputVar = var(outputs,0,2);

R_squar = 1-(errVar/inputVar);

clear inputs outputs netOut err errVar inputVar
%% coparing NN outputs to Matsuoka's estimation:
% attention: check how you defined the NN inputs.
% define the inputs order as: 'tau','b','a','s'
% define the output as: 'frequency'
tau = sampl(1,:);
T = 5*tau;
b = sampl(2,:);
a = sampl(3,:);
cond1 = 1.2 < a; % 1+tau/T < a
cond2 = a < 1+b; % a<1+b

good_ids = find(cond1 & cond2);
tau = tau(good_ids);
T = 5*tau;
b = b(good_ids);
a = a(good_ids);

w_n = (1./T).*sqrt(((tau+T).*b-a.*tau)./(a.*tau));
w_n_fromCode = targ(:,good_ids);
w_n_fromNN = net(sampl(:,good_ids));
errMSE_estVScode = immse(w_n,w_n_fromCode); % calculate MSE error
figure;
subplot(2,1,1);
scatter(w_n_fromCode,w_n);
xlabel('\omega_{n} from code'); ylabel('\omega_{n} from Matsuokas estimation');
subplot(2,1,2);
scatter(w_n_fromCode,w_n); axis([0,8,0,8]);
xlabel('\omega_{n} from code'); ylabel('\omega_{n} from Matsuokas estimation');

figure;
scatter(w_n_fromCode,w_n_fromNN);
xlabel('\omega_{n} from code'); ylabel('\omega_{n} from NN');

figure;
subplot(2,1,1);
scatter(w_n,w_n_fromNN);
xlabel('\omega_{n} from Matsuokas estimation'); ylabel('\omega_{n} from NN');
subplot(2,1,2);
scatter(w_n,w_n_fromNN); axis([0,8,0,8]);
xlabel('\omega_{n} from Matsuokas estimation'); ylabel('\omega_{n} from NN');

clear tau T b a w_n w_n_fromCode w_n_fromNN
clear cond1 cond2
%% checking populations in each "cloud" or "Lobe"
% attention: check how you defined the NN inputs and outputs.
% define the inputs order as: 'tau','b','a','s'
% define the output as: 'frequency'
w_n_fromCode = targ;
w_n_fromNN = net(sampl);

slope = 2.949/1.954; % the slope of the line which seperationg the two lobes (x,y taken by hand);
% make slope '1' to get the 45[deg] line (where NNout=targ)

upper_lobe = find(w_n_fromNN > slope*w_n_fromCode);
lower_lobe = find(w_n_fromNN < slope*w_n_fromCode);

figure; hold on;
scatter(w_n_fromCode(1,upper_lobe),w_n_fromNN(1,upper_lobe),'r');
scatter(w_n_fromCode(1,lower_lobe),w_n_fromNN(1,lower_lobe),'b');
xlabel('\omega_{n} from code'); ylabel('\omega_{n} from NN');

tau = sampl(1,:);
b = sampl(2,:);
a = sampl(3,:);
a_b_ratio = a./b;
figure;
subplot(2,2,1);
histogram(a(1,upper_lobe),100);
xlabel('a'); title('Hist of a in the upper lobe');
subplot(2,2,2);
histogram(a_b_ratio(1,upper_lobe),100);
xlabel('a/b'); title('Hist of a/b in the upper lobe');
subplot(2,2,3);
histogram(a(1,lower_lobe),100);
xlabel('a'); title('Hist of a in the lower lobe');
subplot(2,2,4);
histogram(a_b_ratio(1,lower_lobe),100);
xlabel('a/b'); title('Hist of a/b in the lower lobe');

figure;
subplot(2,1,1);
histogram(tau(1,upper_lobe),100);
xlabel('\tau'); title('Hist of \tau in the upper lobe');
subplot(2,1,2);
histogram(tau(1,lower_lobe),100);
xlabel('\tau'); title('Hist of \tau in the lower lobe');% trying to identify with Matsuoka's condition 
cond1 = 1.2 < a; % 1+tau/T < a
cond2 = a < 1+b; % a<1+b
cond1_ind = find(cond1 & ~cond2);
cond2_ind = find(cond2 & ~cond1);
both_ind = find(cond1 & cond2);
none_ind = find(~cond1 & ~cond2);

figure; hold on
scatter(w_n_fromCode(1,cond1_ind),w_n_fromNN(1,cond1_ind),'b');
scatter(w_n_fromCode(1,cond2_ind),w_n_fromNN(1,cond2_ind),'g');
scatter(w_n_fromCode(1,both_ind),w_n_fromNN(1,both_ind),'k');
scatter(w_n_fromCode(1,none_ind),w_n_fromNN(1,none_ind),'r');
xlabel('\omega_{n} from code'); ylabel('\omega_{n} from NN');
legend('blue: 1+\tau/T<a but no a<1+b','green: a<1+b but no 1+\tau/T<a',...
    'black: both condition exist','red: none of them');

%% plotting one gene from each group:

randUP_gene = randsample(upper_lobe,1); %the index taken from NN data group
randLOW_gene = randsample(lower_lobe,1);%the index taken from NN data group

period_fromCode_UP = 1/targ(1,randUP_gene); % taking NN estimation
period_fromNN_UP = 1/net(sampl(:,randUP_gene));
period_fromCode_LOW = 1/targ(1,randLOW_gene);
period_fromNN_LOW = 1/net(sampl(:,randLOW_gene));

upperGroupCPG = ids(randUP_gene); % convert the index to the indexing of all data ('results')
lowerGroupCPG = ids(randLOW_gene); % convert the index to the indexing of all data ('results')

MML = MatsuokaML(); % calling Matsuoka class
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01; % 0.01
MML.tEnd = 15; % 15

[outUP, ~, signalUP] = MML.runSim(results(upperGroupCPG).seq);
[outLOW, ~, signalLOW] = MML.runSim(results(lowerGroupCPG).seq);


figure; 
subplot(2,2,1);
plot(signalUP.T,signalUP.y);
xlabel('time [sec]'); ylabel('y_i');
title('neurons output from the upper group');
subplot(2,2,2);
plot(signalLOW.T,signalLOW.y);
xlabel('time [sec]'); ylabel('y_i');
title('neurons output from the lower group');
subplot(2,2,3);
plot(signalUP.T,signalUP.signal'); %axis([20,22,-1,1]); 
xlabel('time [sec]'); ylabel('y_2-y_1');
title({'CPG output from the upper group',['geneNum: ',...
    num2str(upperGroupCPG)],...
    [' period_{code} = ',num2str(outUP.periods),'[sec]   '...
    '   period_{NN} ',num2str(period_fromNN_UP),'[sec]']});
subplot(2,2,4);
plot(signalLOW.T,signalLOW.signal');% axis([20,22,-2,2])
xlabel('time [sec]'); ylabel('y_2-y_1');
title({'CPG output from the lower group',['geneNum: ',...
    num2str(lowerGroupCPG)],...
    [' period_{code} = ',num2str(outLOW.periods),'[sec]   '...
    '   period_{NN} ',num2str(period_fromNN_LOW),'[sec]']});

%% rum one gene (from the upper lobe) many times and check distribution of the period
N = 1000;
periodsSim = zeros(1,N);

for i=1:N
    [outUP, ~, ~] = MML.runSim(results(upperGroupCPG).seq);
    periodsSim(1,i) = outUP.periods;
    clear outUP
end

figure;
subplot(2,1,1);
scatter(1:N,periodsSim);
subplot(2,1,2)
histogram(periodsSim,10);
