
clear all; close all; clc;

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
%% Phase 2 - load data to arrays:
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

%% Phase 3 - load data to NN training:
fileName_for_train = 'MatsRandomRes_2-4_02_2017.mat';
myCode = myCode(fileName_for_train,parametersCells,targetCells,seqOrder,4);
myCode.sizeOfCPG = 4;

%% manually:
netFile = uigetfile;
netemp = load(netFile);
net = netemp.NN;   clear netemp netFile

[percent_osc_new,conv_in_range,accuracy] = ...
        NN_GA_perf(sampl,targ,results_old,fileName,seqOrder,caseNum,perORfreq);
%% automated
percent_osc_new = zeros(1,10);
conv_in_range = zeros(1,10);
accuracy = zeros(1,10);
NN_perf = zeros(1,10);

% what NN case do we want:
caseNum = 1:10; % options: '1'-'10'

for i=1:10
    
    disp([' at i = ',num2str(i)]);
    myCode.disp_information = false;
    myCode = myCode.Set('NN',[10],500);
    myCode = myCode.trainNN(0,0);
    net = myCode.NN.net;
    
    NN_perf(1,i) = myCode.NN.MSE_test_perf;
    
    [percent_osc_new(1,i),conv_in_range(1,i),accuracy(1,i)] = ...
        NN_GA_perf(net,sampl,results_old,fileName,...
        seqOrder,caseNum(1,i),perORfreq);
end

%% Plot

figure;
bar(NN_perf);
title('NN with 10 hidden neurons and different inputs outputs');
xlabel('#case');
ylabel('MSE perf');

figure;
bar([percent_osc_new;conv_in_range;accuracy]);
title('NN with 10 hidden neurons and different inputs outputs');
xlabel('#case');
legend('percent osc new','conv in range','accuracy');