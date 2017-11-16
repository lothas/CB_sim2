%% NN results:
clear all; close all; clc

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;

% % % data with many CPG's that oscillates in range:
results_fileName = 'MatsRandomRes_4Neurons_Large_b_Large_W_All_osc';

% % data with small amount of CPGs that oscillate in range:
% results_fileName = 'MatsRandomRes_4Neurons_4Paper.mat';

% % % % data with tau_ratio=12 and (0.2 < b < 2.5):
% % change b_max:
% MML.Gen.Range(2,2) = 2.5; % the class will filter genes that are not in the new range.
% results_fileName = 'MatsRandomRes_4Neurons_4Paper_tau_ratio_equalTo_12_added_b_4Paper1.mat';
% results_fileName = 'MatsRandomRes_4Neurons_4Paper_narrower_W_range';

% % change tau_a/tau_r to 12 (instead of 5)
MML.Sim.Con.tau_ratio = 12;

NNs_4paper = NNs_4paper(results_fileName,MML);

%% plot histograms of  osc_param Vs n-osc_param:
paramName = 'tau';
norm_flag = true;
NNs_4paper.plot_oscParam_vs_NoscParam_hist(paramName,norm_flag)
NNs_4paper.plot_oscParam_vs_oscInRangeParam_hist(paramName,norm_flag);

paramName = 'b';
norm_flag = false;
NNs_4paper.plot_oscParam_vs_NoscParam_hist(paramName,norm_flag);
NNs_4paper.plot_oscParam_vs_oscInRangeParam_hist(paramName,norm_flag);

paramName = 'w_{12}';
norm_flag = false;
NNs_4paper.plot_oscParam_vs_NoscParam_hist(paramName,norm_flag);
NNs_4paper.plot_oscParam_vs_oscInRangeParam_hist(paramName,norm_flag);
%% plot 2D histogram:
norm_flag = true;
figure;
NNs_4paper.plot_2D_hist({'tau','b'},norm_flag)

periods = horzcat(NNs_4paper.results(NNs_4paper.osc_ids).periods);
figure;
histogram(periods,100,'Normalization','pdf');
xlabel('period [sec]');
title('histogram of the CPGs period');

%% rescale tau
NNs_4paper = NNs_4paper.rescale_CPGs();

tau_before = NNs_4paper.seq(:,1);
tau_after = NNs_4paper.tau_rescaled;

% norm results with min,max:
tau_before = NNs_4paper.norm_min_max(tau_before,'tau');
tau_after = NNs_4paper.norm_min_max(tau_after,'tau');

NNs_4paper.hist_compare(tau_before,tau_after,'tau',...
    20,{'tau','tau_{rescaled}'},'plot')

%% check NN on training data:
close all; clc;

caseNum = 7;
% get the names of the training parameters:
[Inputs_names,Targets_names] =...
    NNs_4paper.check_NN_case(caseNum,'period');

% Inputs_names = {'tau','b','w_{12}','w_{13}','w_{14}',...
%         'w_{21}','w_{23}','w_{24}',...
%         'w_{31}','w_{32}','w_{34}',...
%         'w_{41}','w_{42}','w_{43}'};
% Targets_names = {'period'};


architecture = [20];

NNs_4paper = NNs_4paper.train_and_test(Inputs_names,Targets_names,...
    architecture,'NN',1);

% NNs_4paper.train_and_test(Inputs_names,Targets_names,...
%     architecture,'MoE colaboration',1);
% 
% NNs_4paper.train_and_test(Inputs_names,Targets_names,...
%     architecture,'MoE hard',1);
% 
% NNs_4paper.train_and_test(Inputs_names,Targets_names,...
%     architecture,'MoE soft',1);


%% clustering attempt:

close all; clc

n = 3;

inputs = NNs_4paper.Inputs_train;
targets = NNs_4paper.Targets;
x = inputs;

% net_clust = selforgmap([n,n]);
net_clust = competlayer(n*n);
net_clust = train(net_clust,x);
view(net_clust)
y = net_clust(x);
classes = vec2ind(y);

figure;
targets_clustered = cell(1,n*n);
for i=1:n*n 
    subplot(n,n,i)
    ind = find(classes == i);
    histogram(targets(1,ind));
    
    targets_clustered{1,i} = targets(1,ind);
end

figure;
net_reg = cell(1,n*n);
outputs_clustered = cell(1,n*n);
for i=1:n*n
    disp(['at i=',num2str(i)]);
    
    ind = find(classes == i);
    
    tempNet = fitnet(15);
    tempNet.trainFcn = 'trainscg';
    tempNet.trainParam.showWindow = 0; 
    [tempNet, tr] = train(tempNet,...
        inputs(:,ind),...
        targets(:,ind),...
        'useParallel','yes','useGPU','yes');

    net_reg{1,i} = tempNet;
    
    outputs_clustered{1,i} = tempNet(inputs(:,ind));
    RMSE = sqrt(immse(targets(:,ind),outputs_clustered{1,i}));
    
    subplot(n,n,i);
    histogram(targets(:,ind),100); hold on;
    histogram(outputs_clustered{1,i},100);
    title(['RMSE = ',num2str(RMSE)]);
    
    clear tempNet RMSE tr
    
end

figure;
for i=1:n*n
    subplot(n,n,i);
    scatter(targets_clustered{1,i},outputs_clustered{1,i});
    title([num2str(length(targets_clustered{1,i})),' sample']);
end


%% PCA trying:

% data = [targets;inputs];
% 
% categories = {'b','periods','tau',...
%     'w_{12}','w_{13}','w_{14}',...
%     'w_{21}','w_{23}','w_{24}',...
%     'w_{31}','w_{32}','w_{34}',...
%     'w_{41}','w_{42}','w_{43}'};

data = [inputs];

categories = {'periods','tau',...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};

figure;
boxplot(data','orientation','horizontal','labels',categories)

% Check the pairwise correlation between the variables:
C = corr(data',data');

figure;
imagesc(C);

w = 1./var(data');
[wcoeff,score,latent,tsquared,explained] = pca(data',...
'VariableWeights',w);

figure()
plot(score(:,1),score(:,2),'+')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')

figure()
pareto(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')

% Transform the coefficients so that they are orthonormal
coefforth = inv(diag(std(data')))*wcoeff;

figure;
whichCompo = 1:2;
ind = randsample(size(data,2),1000);
biplot(coefforth(:,whichCompo),...
    'scores',score(ind,whichCompo),...
    'varlabels',categories);
%% train NN 5 times and collect statistics about the Perf:
close all; clc;
caseNum = 7;
% get the names of the training parameters:
[Inputs_names,Targets_names] =...
    NNs_4paper.check_NN_case(caseNum,'period');

% Inputs_names = {'tau','b','w_{12}','w_{13}','w_{14}',...
%         'w_{21}','w_{23}','w_{24}',...
%         'w_{31}','w_{32}','w_{34}',...
%         'w_{41}','w_{42}','w_{43}'};
% Targets_names = {'period'};

architecture = {[20,20]};
numOfRepeats = 5;

train_RMSE = zeros(numOfRepeats,length(architecture));
valid_RMSE = zeros(numOfRepeats,length(architecture));
test_RMSE = zeros(numOfRepeats,length(architecture));

train_R2 = zeros(numOfRepeats,length(architecture));
valid_R2 = zeros(numOfRepeats,length(architecture));
test_R2 = zeros(numOfRepeats,length(architecture));

train_slope = zeros(numOfRepeats,length(architecture));
valid_slope = zeros(numOfRepeats,length(architecture));
test_slope = zeros(numOfRepeats,length(architecture));

for i=1:length(architecture)
    disp(['at NN arch :  ',num2str(architecture{1,i})]);
    for j=1:5
        disp(['    trail num #',num2str(j)])
        NNs_4paper = NNs_4paper.train_and_test(Inputs_names,Targets_names,...
        architecture{1,i},'NN',0);

        train_RMSE(j,i) = NNs_4paper.NN.train_RMSE;
        valid_RMSE(j,i) = NNs_4paper.NN.valid_RMSE;
        test_RMSE(j,i) = NNs_4paper.NN.test_RMSE;

        train_R2(j,i) = NNs_4paper.NN.train_R2;
        valid_R2(j,i) = NNs_4paper.NN.valid_R2;
        test_R2(j,i) = NNs_4paper.NN.test_R2;

        train_slope(j,i) = NNs_4paper.NN.train_slope;
        valid_slope(j,i) = NNs_4paper.NN.valid_slope;
        test_slope(j,i) = NNs_4paper.NN.test_slope;
    end
end

train_RMSE_mean = mean(train_RMSE,1);
valid_RMSE_mean = mean(valid_RMSE,1);
test_RMSE_mean = mean(test_RMSE,1)

train_R2_mean = mean(train_R2,1);
valid_R2_mean = mean(valid_R2,1);
test_R2_mean = mean(test_R2,1)

train_slope_mean = mean(train_slope,1);
valid_slope_mean = mean(valid_slope,1);
test_slope_mean = mean(test_slope,1)

% figure;
% boxplot(train_RMSE,'Colors',[0,0,128]./256); hold on;
% boxplot(valid_RMSE,'Colors',[34,139,34]./256);
% boxplot(test_RMSE,'Colors',[178,34,34]./256);
% xlabel('NN arcith');
% ylabel('RMSE');
% grid minor
% title('RMSE over hidden neurons num');
% 
% figure;
% boxplot(train_R2,'Colors',[0,0,128]./256); hold on;
% boxplot(valid_R2,'Colors',[34,139,34]./256);
% boxplot(test_R2,'Colors',[178,34,34]./256);
% xlabel('NN arcith');
% ylabel('R^2');
% grid minor
% title('R^2 over hidden neurons num');
% 
% figure;
% boxplot(train_slope,'Colors',[0,0,128]./256); hold on;
% boxplot(valid_slope,'Colors',[34,139,34]./256);
% boxplot(test_slope,'Colors',[178,34,34]./256);
% xlabel('NN arcith');
% ylabel('reggresion graph slope');
% grid minor
% title('slope over hidden neurons num');

% figure;hold on;
% H(1) = shadedErrorBar(x, y, {@mean, @(x) 2*std(x)  }, '-r', 0);
% H(2) = shadedErrorBar(x, y, {@mean, @(x) 1*std(x)  }, '-m', 0);
% H(3) = shadedErrorBar(x, y, {@mean, @(x) 0.5*std(x)}, {'-b', 'LineWidth', 2}, 0);
% 
% legend([H(3).mainLine, H.patch], ...
%     '\mu', '2\sigma', '\sigma', '0.5\sigma', ...
%     'Location', 'Northwest');
%% save data to CSV file:
caseNum = 9;
fileName = 'data_to_CSV_case_7_tauRatio_12_';
NNs_4paper.save_data_CSV(caseNum,fileName);
%% correlation between targets and NN outputs:
% a = [tau_out_on_train_set;b_out_on_train_set];
% b = [tau_training;b_training];
% [RHO,PVAL] = corr(a',b')