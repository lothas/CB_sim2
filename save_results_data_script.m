
close all; clc; clear all;

% % Load peridos oscillating and periods which oscillates in range:
load('MatsRandomRes_4Neurons_Large_b_Large_W_All_osc1.mat','results','header');
results1 = results;    clear results
load('MatsRandomRes_4Neurons_Large_b_Large_W_All_osc2.mat','results');
results2 = results;    clear results
results = [results1,results2];  clear results1 results2


% header = sprintf('tau ratio is equal to 12 \n');
% header = [header,sprintf('data is for 2N symmetric CPG case \n')];
% header = [header,sprintf('seq Order: \n')];
% header = [header,sprintf('"tau" ,"b", "c", NR, "a", NR... \n')];
% header = [header,sprintf('"NR" - not relevnt param \n')];
% header = [header,sprintf('b in range (0.2,2.5) \n')];
% header = [header,sprintf('All oscillating in range')];

save('MatsRandomRes_4Neurons_Large_b_Large_W_All_osc.mat',...
    'header','results');

%% save the "good" samples:

header = sprintf('tau ratio is equal to 12 \n');
header = [header,sprintf('data is for 4N general CPG case \n')];
header = [header,sprintf('seq Order: \n')];
header = [header,sprintf('"tau","b","c_1","c_2","c_3","c_4" \n')];
header = [header,sprintf('"w_{12}","w_{13}","w_{14}","w_{21}","w_{23}","w_{24}" \n')];
header = [header,sprintf('"w_{31}","w_{32}","w_{34}","w_{41}","w_{42}","w_{43}" \n')];
header = [header,sprintf('b in range (0.2,10) \n')];
header = [header,sprintf('W_ij in range (0,10) \n')];

% results = NNs_4paper.results(NNs_4paper.osc_ids);
results = NNs_4paper.results(NNs_4paper.osc_ids);

save('MatsRandomRes_4Neurons_Large_b_Large_W_All_osc2.mat',...
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