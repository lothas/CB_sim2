
clear all; close all; clc;

%NOTE: if you want to run automatically on all cases Run Phase 2!!!
%% Create Genome File
genome_file = 'MatsuokaGenome.mat';
nAnkle = 1; % Number of ankle torques
nHip = 1;   % Number of hip torques
maxAnkle = 10;%20;   % Max ankle torque
maxHip = 10;%20;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;
Mw = 10*ones(1,12);
mw = 0.1*ones(1,12);

%     %%%%%%%%%%%% For the 4-neuron case!!!
%     % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
Keys = {'\tau_r', 'beta', 'amp',   'weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
              1 ,      1,     4,          12,        1 ,       2*N ,            0 };
Range = {  0.02 ,    0.2,  mamp,          mw,   -0.001 ,  -0.2*Mamp; % Min
           0.25 ,   10.0,  Mamp,          Mw,    0.001 ,   0.2*Mamp}; % Max

MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

clear all;
%% Initialize machine learning object for Matsuoka analysis
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01; % 0.01
MML.tEnd = 30;
nPlotSamples = 0; % 10;
% Turn off findpeaks warning
warning('off','signal:findpeaks:largeMinPeakHeight');

%% Phase 0 - Run 2000 Matsuoka simulations with different parameters
filename1 = 'MatsRandomRes_Test_NN_allSims.mat';
nSamples = 5000; % run 2000, take 500 with period
MML.runRandomSims(nSamples, filename1);
clear nSamples

%% Phase 1- load 600 random Non-oscillating CPG's
load('MatsRandomRes_Test_NN_allSims.mat');
clear periods filename1 id_conv id_per nPlotSamples nSims

periods = horzcat(results(:).periods);
periods = periods(1,:);

% osc_ids = find(~isnan(periods)); % find oscilatorry
osc_ids = find(isnan(periods)); % find NOT oscilatorry
cpg_toUse = randsample(osc_ids,600);
results1 = results;             clear results;
results = results1(cpg_toUse);  clear results1;

save('MatsRandomRes_Test_NN_Sims2use.mat','results')
clear periods osc_ids results

%% Phase 2: automated
clear all; clc
% what NN case do we want:
caseNum = [1,3,5,7,9]; % options: '1'-'10'
% caseNum = [1,7]; % options: '1'-'10'
N = length(caseNum);

% NN information:
hidN_num = 20;

% number of repeats:
repeat_num = 5;

seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};
fileName_4GA = 'MatsRandomRes_Test_NN_allSims.mat';
% fileName_for_train = 'MatsRandomRes_2-4_02_2017.mat';
fileName_for_train = 'MatsRandomRes_all_from_1-2_2017.mat';


for i=1:N
    
    disp([' testing NN_case #',num2str(caseNum(i))]);
    
    NN_perf_tau = nan(1,repeat_num);
    NN_perf_b = nan(1,repeat_num);
    percent_osc_new = zeros(1,repeat_num);
    conv_in_range = zeros(1,repeat_num);
    accuracy = zeros(1,repeat_num);
    
    % make input and outputs names according to what NN case do we want:
    %   Testing case:
    [parametersCells_test,targetCells_test] = ...
        check_NN_case_for_paper(caseNum(1,i),'period_desired');
    
    % Load data for GA testing (500-600 samples):
    myCode_test = myCode(fileName_4GA,parametersCells_test,...
        targetCells_test,seqOrder,4);
    myCode_test.sizeOfCPG = 4;
%     [results_old,sampl,~] = myCode.NN_prepare_badCPGs_to_NN(600);
    [results_old,sampl,~] = myCode_test.NN_prepare_badCPGs_to_NN(100);
    
    % make input and outputs names according to what NN case do we want:
    %   Training case:
    [parametersCells_train,targetCells_train] = ...
        check_NN_case_for_paper(caseNum(1,i),'period');
    
    % Load data for NN training (lots of samples):
    myCode_train = myCode(fileName_for_train,parametersCells_train,...
        targetCells_train,seqOrder,4);
    myCode_train.sizeOfCPG = 4;
    myCode_train.train_ratio = 0.7;
    myCode_train.valid_ratio = 1-0.1-0.7; % '0.1'=test group size
    
    for j=1:repeat_num % repeat 5 times for statistics reasons
        disp(['    iter #',num2str(j)]);
        
        % Train NN:
        myCode_train.disp_information = false;
        myCode_train = myCode_train.Set('NN',[hidN_num],1000);
        myCode_train = myCode_train.trainNN(0,0);
        net = myCode_train.NN.net;

        switch caseNum(1,i)
            case {1,2,3,4}
                NN_perf_tau(1,j) = myCode_train.NN.MSE_test_perf;
            case {5,6,7,8}
                NN_perf_b(1,j) = myCode_train.NN.MSE_test_perf;
            case {9,10}
                test_targ = myCode_train.targ_test;
                test_NNout = net(myCode_train.sampl_test);
                % calculate MSE error on 'tau':
                [NN_perf_tau(1,j),~] = ...
                    myCode_train.NN_perf_calc(test_targ(1,:),test_NNout(1,:),...
                    0,0,'test' );
                % calculate MSE error on 'b':
                [NN_perf_b(1,j),~] = ...
                    myCode_train.NN_perf_calc(test_targ(2,:),test_NNout(2,:),...
                    0,0,'test' );
        end



        [percent_osc_new(1,j),conv_in_range(1,j),accuracy(1,j)] = ...
            NN_GA_perf(net,sampl,results_old,seqOrder,caseNum(1,i));

        % Shuffle data to get new training group:
        myCode_train = myCode_train.shuffle_samples('completeShuffle');
    end
    
    NNs_cases_perf(i).case = caseNum(1,i);
    NNs_cases_perf(i).NN_perf_tau = NN_perf_tau;
    NNs_cases_perf(i).NN_perf_b = NN_perf_b;
    NNs_cases_perf(i).percent_osc_new = percent_osc_new;
    NNs_cases_perf(i).conv_in_range = conv_in_range;
    NNs_cases_perf(i).accuracy = accuracy;

end

fileName = 'testing_NNs_cases.mat';
save(fileName,'NNs_cases_perf');
%% Phase 4- plot and save results:
fileName = 'testing_NNs_cases.mat';
load(fileName,'NNs_cases_perf');

hidN_num = 20;

caseNum = (vertcat(NNs_cases_perf(:).case))';
NN_perf_tau = vertcat(NNs_cases_perf(:).NN_perf_tau);
NN_perf_b = vertcat(NNs_cases_perf(:).NN_perf_b);
percent_osc_new = vertcat(NNs_cases_perf(:).percent_osc_new);
conv_in_range = vertcat(NNs_cases_perf(:).conv_in_range);
accuracy = vertcat(NNs_cases_perf(:).accuracy);

% calc mean and stdev
NN_perf_tau_mean = (mean(NN_perf_tau,2))';
NN_perf_b_mean = (mean(NN_perf_b,2))';
percent_osc_new_mean = (mean(percent_osc_new,2))';
conv_in_range_mean = (mean(conv_in_range,2))';
accuracy_mean = (mean(accuracy,2))';

% calc mean and stdev
NN_perf_tau_std = (std(NN_perf_tau,[],2))';
NN_perf_b_std = (std(NN_perf_b,[],2))';
percent_osc_new_std = (std(percent_osc_new,[],2))';
conv_in_range_std = (std(conv_in_range,[],2))';
accuracy_std = (std(accuracy,[],2))';

% plot
figure; hold on
h1 = bar(caseNum,NN_perf_tau_mean);
h1.FaceAlpha = 0.5;
errorbar(caseNum,NN_perf_tau_mean,NN_perf_tau_std,'.')
title(['NN with ',num2str(hidN_num),...
    ' hidden neurons and different inputs outputs']);
xlabel('#case');
ylabel('MSE perf on \tau');
grid minor

figure; hold on
bar(caseNum,NN_perf_b_mean);
errorbar(caseNum,NN_perf_b_mean,NN_perf_b_std,'.')
title(['NN with ',num2str(hidN_num),...
    ' hidden neurons and different inputs outputs']);
xlabel('#case');
ylabel('MSE perf on b');
grid minor; hold off;

% get class "myCode" for use of funtion "bars..."
[parametersCells,targetCells] = check_NN_case_for_paper(1,'period_desired');
fileName_4GA = 'MatsRandomRes_Test_NN_allSims.mat';
myCode = myCode(fileName_4GA,parametersCells,targetCells,seqOrder,4);
myCode.sizeOfCPG = 4;

means = [percent_osc_new_mean;...
    conv_in_range_mean;...
    accuracy_mean];
stdevs = [percent_osc_new_std;...
    conv_in_range_std;...
    accuracy_std];
Names = {' ',' ','case#1',' ',' ',' ',' ','case#3',' '};
label_Y = '';
graph_title = ['NN with ',num2str(hidN_num),...
    ' hidden neurons and different inputs outputs'];
graph_legend = {'converged','conv in range','accuracy'};
myCode.plot_MSE_perf_with_stdev(means,stdevs,...
    Names,label_Y,graph_title,graph_legend)

% % plot
% errBarLocs = caseNum'*ones(1,3) + ones(5,1)*[-0.5, 0, 0.5 ];
% figure; hold on
% h1 = bar(caseNum'*ones(1,3),means');
% h1.FaceAlpha = 0.5;
% errorbar(errBarLocs,means',stdevs','.')
% title(['NN with ',num2str(hidN_num),...
%     ' hidden neurons and different inputs outputs']);
% xlabel('#case');
% ylabel('MSE perf on \tau');
% grid minor
