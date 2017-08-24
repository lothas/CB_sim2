
close all; clc; clear all;

% % Load peridos oscillating and periods which oscillates in range:
load('MatsRandomRes_4Neurons_4Paper_tau_ratio_equalTo_12_more_natural_A','results','MML');
results1 = results;    clear results
load('MatsRandomRes_4Neurons_4Paper_tau_ratio_equalTo_12_more_natural_B','results');
results2 = results;    clear results
load('MatsRandomRes_4Neurons_4Paper_tau_ratio_equalTo_12_more_natural_C','results');
results3 = results;    clear results
results = [results1,results2,results3];  clear results1 results2 results3


header = sprintf('tau ratio is equal to 12 \n');
header = [header,sprintf('b in range (0.2,2.5)')];

save('MatsRandomRes_4Neurons_4Paper_tau_ratio_equalTo_12_added_b_All_samples.mat',...
    'header','results','MML');

%% save the "good" samples:

header = sprintf('tau ratio is equal to 12 \n');
header = [header,sprintf('sample with b in range (0.2,2.5) \n')];
header = [header,sprintf('filtered the non-osc samples')];

results = NNs_4paper.results(NNs_4paper.osc_ids);

save('MatsRandomRes_4Neurons_4Paper_tau_ratio_equalTo_12_more_osc_samples.mat',...
    'header','results');

%% combaining the filtered CPGs
close all; clc; clear all;

% % Load peridos oscillating and periods which oscillates in range:
load('MatsRandomRes_4Neurons_4Paper_tau_ratio_equalTo_12_added_b_4Paper.mat','results');
results1 = results;    clear results
load('MatsRandomRes_4Neurons_4Paper_tau_ratio_equalTo_12_more_natural_H.mat','results');
results2 = results;    clear results

results = [results1,results2];  
clear results1 results2


header = sprintf('tau ratio is equal to 12 \n');
header = [header,sprintf('b in range (0.2,2.5) \n')];
header = [header,sprintf('filtered the non-osc samples and then added around 40k non-osc samples for later \n')];
header = [header,sprintf('originally, the ratio between osc to non-osc was 1 (osc) to 6 (n-osc')];

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;
MML.Sim.Con.tau_ratio = 12;
MML.Gen.Range(2,2) = 2.5; % the class will filter genes that are not in the new range.


save('MatsRandomRes_4Neurons_4Paper_tau_ratio_equalTo_12_added_b_4Paper1.mat',...
    'header','results','MML');

%% just add header to file:
header = sprintf('tau ratio is equal to 2 \n');
header = [header,sprintf('half of the samples are oscillatory \n')];
header = [header,sprintf('and the other half where rescaled to the desired period range \n')];
header = [header,sprintf('originally, the ratio between osc to non-osc was 1 (osc) to 2 (n-osc)')];

save('MatsRandomRes_4Neurons_4Paper_for_MOOGA_try.mat',...
    'header','results','MML');