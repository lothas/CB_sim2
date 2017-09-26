
%% NN testing Script:
% IMPORTANT: dont forget to load the right genome file and to uptade
% 'MatsuokaML.m' to the right rettings
% 

clear all; close all; clc;

%%

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau' ,'b', 'c', 'NR', 'a',...
    'k_tau','k_{c}'};
% "NR" - not relevnt param 

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;

% % file name for uploading:
% results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_only_osc_1.mat'};

results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_only_osc_1.mat',...
    'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_only_osc_3.mat',...
    'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_only_osc_4.mat'};

%% Load data:
results = load_results(results_fileName);


%% get and filter periods:
% % get oscillating:
[results,periods,seq,~] = get_CPGs(results,'osc','2N',MML);

% % % get oscillating in period range:
% [results,periods,seq,~] = get_CPGs(results,'osc_in_per_range','2N',MML);


%% Neural Network #1:
input_names_net1 = {'b','tau','a'};
output_names_net1 = {'periods'};
[sampl_net1,targ_net1] = ...
    prepare_NN_data(input_names_net1,output_names_net1,...
    seqOrder,seq,periods);

architecture = [10];

% minTar = min(targ);
% maxTar = max(targ);
% targ_norm = (targ - minTar) ./ (maxTar - minTar);
targ_norm_net1 = targ_net1; % remove if not necessary

net1 = fitnet(architecture);
% net1 = feedforwardnet(architecture);
net1.trainFcn = 'trainbr';

% % bound the NN output with 'logsig' output function:
net1.layers{2,1}.transferFcn = 'logsig';
net1.output.processParams{1,2}.ymin = 0;

net1.divideParam.trainRatio = 0.5;
net1.divideParam.valRatio = 0.35;
net1.divideParam.testRatio = 0.15;

net1.trainParam.showWindow = 1; 
net1.trainParam.showCommandLine = 1;
net1.trainParam.epochs = 1000;
% net1.performParam.normalization = 'percent'; %It must be 'none', 'standard' or 'percent'

[net1, tr1] = train(net1, sampl_net1, targ_norm_net1);

% NN_out = minTar + (net1(sampl) * (maxTar - minTar));

clear NN_out

%% Neural Network #2:
input_names_net2 = {'tau','a','periods'};
output_names_net2 = {'b'};
[sampl_net2,targ_net2] = ...
    prepare_NN_data(input_names_net2,output_names_net2,...
    seqOrder,seq,periods);

architecture = [10];

net2 = fitnet(architecture);
% net2 = feedforwardnet(architecture);
net2.trainFcn = 'trainbr';

net2.layers{2,1}.transferFcn = 'logsig';
net2.output.processParams{1,2}.ymin = 0;

net2.divideParam.trainRatio = 0.5;
net2.divideParam.valRatio = 0.35;
net2.divideParam.testRatio = 0.15;

net2.trainParam.showWindow = 1; 
net2.trainParam.showCommandLine = 1;
net2.trainParam.epochs = 1000;
% net2.performParam.normalization = 'percent'; %It must be 'none', 'standard' or 'percent'

[net2, tr2] = train(net2, sampl_net2, targ_net2);

clear NN_out

%%
close all
% % % % % 2) Test with external test data (of osc in range CPGs)
% results_test = load_results({'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_osc_inRange_test_group1.mat'});
% [~,~,seq_test,~] = get_CPGs(results_test,'osc_in_per_range','2N',MML);

% % % % % 3) Test with external test data (of osc CPGs)
% results_test = load_results({'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_only_osc_test_group1.mat'});
% [~,~,seq_test,~] = get_CPGs(results_test,'osc','2N',MML);

% % % 4) test on rande seq #1:
seq_test = (MML.Gen.RandSeq(1000))';

% get NN inputs:
[NN1_in_test,~] = ...
    prepare_NN_data(input_names_net1,output_names_net1,...
    seqOrder,seq_test,0);
NN1_out_test = net1(NN1_in_test);
    
% get NN inputs:
[NN2_in_test,targ_test] = ...
    prepare_NN_data(input_names_net2,output_names_net2,...
    seqOrder,seq_test,NN1_out_test);
NN2_out_test = net2(NN2_in_test);
% % % get NN output:
% NN_out_test = minTar + (net(NN_in) * (maxTar - minTar));

plotregression(targ_test,NN2_out_test)

figure;
histogram(NN1_out_test,'Normalization','pdf');
title('histogram of NN1_{output}');
xlabel('NN1_{output}');

figure;
histogram(NN2_out_test,'Normalization','pdf');
title('histogram of NN2_{output}');
xlabel('NN2_{output}');

clear NN_in NN_out desPeriod results_test seq_test