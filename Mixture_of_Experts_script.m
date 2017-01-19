%% Mixture of Experts
clear all; close all; clc;

%% Load data:
load('MatsRandomRes_16_12_2016.mat','results','periods')

%%
ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);
% For 4 neurons:
parametersCells = {'tau','b',...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};

targetCells = {'freq'};
[ sampl,targ ] = prepareData( results,periods,ids,parametersCells,targetCells,0 );

clear results periods ids_period ids_error ids
%% norm inputs
normParams = zeros(size(sampl, 1), 2);
for i = 1:size(sampl, 1)
    feat = sampl(i, :);
    normParams(i, :) = [mean(feat), std(feat)];
    sampl(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
    
end
clear i feat normParams
%% Prepare the NN (="Experts")

% learningRate = 10^-3;
% decay = 0.34;
learningRate = 0.01;
decay = 0.9;
expertCount = 3;
numOfInputs = 14;

% initialize parameters, use first points for m
m = rand(expertCount, numOfInputs);

% define experts:
HiddenN = 10;
net1 = feedforwardnet(HiddenN);
net2 = feedforwardnet(HiddenN);
net3 = feedforwardnet(HiddenN);

net1.trainParam.showWindow = 0; % dont show training window
net2.trainParam.showWindow = 0;
net3.trainParam.showWindow = 0;

maxEphocs = 50;
net1.trainParam.epochs = maxEphocs;
net2.trainParam.epochs = maxEphocs;
net3.trainParam.epochs = maxEphocs;

% initial training (to initilazied the NN)
[net1, tr1] = train(net1, sampl(:,1:50), targ(:,1:50));
[net2, tr2] = train(net2, sampl(:,50:100), targ(:,50:100));
[net3, tr3] = train(net3, sampl(:,100:150), targ(:,100:150));

%% train the experts
numOfIteretions = 20;
expert1GroupSize = zeros(1,numOfIteretions);
expert2GroupSize = zeros(1,numOfIteretions);
expert3GroupSize = zeros(1,numOfIteretions);

for i=1:numOfIteretions
    outMat = [net1(sampl);net2(sampl);net3(sampl)];
    
    errMat = outMat-[targ;targ;targ];
    seMat = errMat.^2; % squar error
    
    [~,best_expert_ind] = min(seMat,[],1);
    
    % clustering to different experts
    cluster1_ind = find(best_expert_ind == 1);
    cluster2_ind = find(best_expert_ind == 2);
    cluster3_ind = find(best_expert_ind == 3);
    
    % training the experts;
    [net1,tr1] = train(net1, sampl(:,cluster1_ind), targ(:,cluster1_ind));
    [net2,tr2] = train(net2, sampl(:,cluster2_ind), targ(:,cluster2_ind));
    [net3,tr3] = train(net3, sampl(:,cluster3_ind), targ(:,cluster3_ind));
    
    expert1GroupSize(1,i) = length(cluster1_ind);
    expert2GroupSize(1,i) = length(cluster2_ind);
    expert3GroupSize(1,i) = length(cluster3_ind);
    
    x_1 = (ones(1,expert1GroupSize(1,i))'*m(1,:));
    x_2 = (ones(1,expert2GroupSize(1,i))'*m(2,:));
    x_3 = (ones(1,expert3GroupSize(1,i))'*m(3,:));
    
%     dm1 = learningRate*sum(sampl(:,cluster1_ind)-x_1',2)/sum(sampl(:,cluster1_ind),2);
%     dm2 = learningRate*sum(sampl(:,cluster2_ind)-x_2',2)/sum(sampl(:,cluster1_ind),2); 
%     dm3 = learningRate*sum(sampl(:,cluster3_ind)-x_3',2)/sum(sampl(:,cluster1_ind),2);
    
    dm1 = learningRate*sum(sampl(:,cluster1_ind)-x_1',2)./sum(sampl(:,cluster1_ind),2);
    dm2 = learningRate*sum(sampl(:,cluster2_ind)-x_2',2)./sum(sampl(:,cluster2_ind),2); 
    dm3 = learningRate*sum(sampl(:,cluster3_ind)-x_3',2)./sum(sampl(:,cluster3_ind),2);
    
    dm=[dm1';dm2';dm3'];
    
    m = m + dm;

    learningRate = learningRate * decay;
end

figure
groupInd = cluster1_ind;
Outputs = net1(sampl);
trOut = Outputs(:,groupInd);
trTarg = targ(groupInd);
plotregression(trTarg,trOut,'Train');

figure
groupInd = cluster2_ind;
Outputs = net2(sampl);
trOut = Outputs(:,groupInd);
trTarg = targ(groupInd);
plotregression(trTarg,trOut,'Train');

figure
groupInd = cluster1_ind;
Outputs = net1(sampl);
trOut = Outputs(:,groupInd);
trTarg = targ(groupInd);
plotregression(trTarg,trOut,'Train');

figure;
plot(1:numOfIteretions,expert1GroupSize,'b'); hold on;
plot(1:numOfIteretions,expert2GroupSize,'r');
plot(1:numOfIteretions,expert3GroupSize,'g');
xlabel('#iteretion');   ylabel('group size [#points]');
legend('1st expert','2nd expert','3th expert');
clear dm1 dm2 dm3 x_1 x_2 x_3

%% clustering the data: devide for experts:
close all;
if false
    % NN classifier for data gating:
    expert1_targs = (best_expert_ind == 1);
    expert2_targs = (best_expert_ind == 2);
    expert3_targs = (best_expert_ind == 3);
    percep_targ = [expert1_targs;expert2_targs;expert3_targs];
    gateNet = patternnet(10); % view(gateNet);
    gateNet = train(gateNet,sampl,percep_targ);
    inx = gateNet(sampl);
end
if true
    % NN classifier for data gating gating using 'g':
    N = size(sampl,2);
    g = zeros(expertCount,N);
    for i=1:N
        input = sampl(:,i);
        g(:,i) = exp(m*input)/sum(exp(m*input));
        g(:,i) = round(g(:,i),4);
    end
    [~,inx] = max(g,[],1);
    clear N
end

bestExpertsInx = [inx;-best_expert_ind];
disp(['the number of train points correctly classify: ',...
    num2str(length(find(~(sum(bestExpertsInx,1))))),...
    ' out of ',num2str(length(bestExpertsInx))]);
%% check train results:

cluster1_train_ind = find(inx == 1);
cluster2_train_ind = find(inx == 2);
cluster3_train_ind = find(inx == 3);

out1 = net1(sampl(:,cluster1_train_ind));
out2 = net2(sampl(:,cluster2_train_ind));
out3 = net3(sampl(:,cluster3_train_ind));

targ1 = targ(:,cluster1_train_ind);
targ2 = targ(:,cluster2_train_ind);
targ3 = targ(:,cluster3_train_ind);

targ_train = [targ1,targ2,targ3];
outM = [out1,out2,out3];

figure;
scatter(targ_train,outM,'b','o'); hold on
% scatter(targ,outMat(1,:),'r','d');
% plot([0, 1],[0, 1],'--k');
xlabel('Targets');  ylabel('NN out');
title('NN regression performance');
plotregression(targ_train,outM,'Train');

plotregression(targ1,out1,'Train');
%% Load test samples:
load('MatsRandomRes_16_12_2016.mat','results','periods')

ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);

targetCells = {'freq'};
[ sampl_test,targ_test ] = prepareData( results,periods,ids,parametersCells,targetCells,0 );

clear results periods ids_period ids_error ids

normParams = zeros(size(sampl_test, 1), 2);
for i = 1:size(sampl_test, 1)
    feat = sampl_test(i, :);
    normParams(i, :) = [mean(feat), std(feat)];
    sampl_test(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
    
end
clear i feat normParams

%% check test results:
N = size(sampl_test,2);
g = zeros(expertCount,N);
outM_test = zeros(1,N);
for i=1:N
    input = sampl_test(:,i);
    g(:,i) = exp(m*input)/sum(exp(m*input));

end

% errMat = outM_test-targ_test;

[~,inx] = max(g,[],1);
cluster1_test_ind = find(inx == 1);
cluster2_test_ind = find(inx == 2);
cluster3_test_ind = find(inx == 3);
out1 =net1(sampl_test(:,cluster1_test_ind));
out2 =net2(sampl_test(:,cluster2_test_ind));
out3 =net3(sampl_test(:,cluster3_test_ind));
targ1_test = targ(:,cluster1_test_ind);
targ2_test = targ(:,cluster2_test_ind);
targ3_test = targ(:,cluster3_test_ind);
targtest = [targ1_test,targ2_test,targ3_test];
outM_test = [out1,out2,out3];
errMat = outM_test-targtest;

figure;
scatter(targ_test,outM_test,'b','o'); hold on
% plot([0, 1],[0, 1],'--k');
xlabel('Targets');  ylabel('NN out');
title('NN regression performance');