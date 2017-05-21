
clear all; close all; clc;

%NOTE: if you want to run automatically on all cases Run Phase 3!!!
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

%% Phase 1- load 800 random oscillating CPG's
load('MatsRandomRes_Test_NN_allSims.mat');
clear periods filename1 id_conv id_per nPlotSamples nSims

periods = horzcat(results(:).periods);
periods = periods(1,:);

osc_ids = find(~isnan(periods));
cpg_toUse = randsample(osc_ids,800);
results1 = results;             clear results;
results = results1(cpg_toUse);  clear results1;

save('MatsRandomRes_Test_NN_Sims2use.mat','results')
clear periods osc_ids results
%% Phase 2 - load data to arrays (RUN manually):
% % %File to load:
fileName = 'MatsRandomRes_Test_NN_Sims2use.mat';
seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};


% what NN case do we want:
caseNum = 3; % options: '1'-'10'
perORfreq = 'freq';
[parametersCells,targetCells] = check_NN_case_for_paper(caseNum,perORfreq);

% % % Load data 
myCode = myCode(fileName,parametersCells,targetCells,seqOrder,4);
myCode.sizeOfCPG = 4;

ind_ids = [myCode.train_ind;myCode.valid_ind;myCode.test_ind];
results_old = myCode.sim_results(ind_ids);
sampl = horzcat(myCode.sampl_train,myCode.sampl_valid,myCode.sampl_test);
targ = horzcat(myCode.targ_train,myCode.targ_valid,myCode.targ_test);

clear myCode caseNum

fileName_for_train = 'MatsRandomRes_2-4_02_2017.mat';
myCode = myCode(fileName_for_train,parametersCells,targetCells,seqOrder,4);
myCode.sizeOfCPG = 4;

netFile = uigetfile;
netemp = load(netFile);
net = netemp.NN;   clear netemp netFile

[percent_osc_new,conv_in_range,accuracy] = ...
        NN_GA_perf(sampl,targ,results_old,fileName,seqOrder,caseNum,perORfreq);
%% Phase 3: automated

% what NN case do we want:
caseNum = 1:10; % options: '1'-'10'
% TODO: fix i=9,10
N = length(caseNum);
percent_osc_new = zeros(1,N);
conv_in_range = zeros(1,N);
accuracy = zeros(1,N);
NN_perf_tau = nan(1,N);
NN_perf_b = nan(1,N);

seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};
fileName_4GA = 'MatsRandomRes_Test_NN_Sims2use.mat';
fileName_for_train = 'MatsRandomRes_2-4_02_2017.mat';


for i=1:N
    
    disp([' at i = ',num2str(i)]);
    
    % make input and outputs names according to what NN case do we want:
%     perORfreq = 'freq';
    perORfreq = 'period';
    [parametersCells,targetCells] = ...
        check_NN_case_for_paper(caseNum(1,i),perORfreq);
    
    % Load data for GA testing (500-600 samples):
    myCode = myCode(fileName_4GA,parametersCells,targetCells,seqOrder,4);
    myCode.sizeOfCPG = 4;
    ind_ids = [myCode.train_ind;myCode.valid_ind;myCode.test_ind];
    results_old = myCode.sim_results(ind_ids);
    sampl = horzcat(myCode.sampl_train,myCode.sampl_valid,myCode.sampl_test);
    targ = horzcat(myCode.targ_train,myCode.targ_valid,myCode.targ_test);
    clear myCode
    
    % Load data for NN training (lots of samples):
    myCode = myCode(fileName_for_train,parametersCells,targetCells,seqOrder,4);
    myCode.sizeOfCPG = 4;
    
    % Train NN:
    myCode.disp_information = false;
    myCode = myCode.Set('NN',[10],1000);
    myCode = myCode.trainNN(0,0);
    net = myCode.NN.net;
    
    switch i
        case {1,2,3,4}
            NN_perf_tau(1,i) = myCode.NN.MSE_test_perf;
        case {5,6,7,8}
            NN_perf_b(1,i) = myCode.NN.MSE_test_perf;
        case {9,10}
            test_targ = myCode.targ_test;
            test_NNout = net(myCode.sampl_test);
            % calculate MSE error on 'tau':
            [NN_perf_tau(1,i),~] = ...
                myCode.NN_perf_calc(test_targ(1,:),test_NNout(1,:),...
                0,0,'test' );
            % calculate MSE error on 'b':
            [NN_perf_b(1,i),~] = ...
                myCode.NN_perf_calc(test_targ(2,:),test_NNout(2,:),...
                0,0,'test' );
    end
    
    
    
    [percent_osc_new(1,i),conv_in_range(1,i),accuracy(1,i)] = ...
        NN_GA_perf(net,sampl,results_old,seqOrder,caseNum(1,i));
    
    clear parametersCells targetCells myCode net results_old sampl targ
end

%% Plot

figure;
bar(NN_perf);
title('NN with 10 hidden neurons and different inputs outputs');
xlabel('#case');
ylabel('MSE perf');
grid minor

figure;
bar([percent_osc_new;conv_in_range;accuracy]');
title('NN with 10 hidden neurons and different inputs outputs');
xlabel('#case');
legend('converged','conv in range','accuracy');
grid minor