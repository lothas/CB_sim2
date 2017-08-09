
clear all; close all; clc

% % Load peridos oscillating and periods which oscillates in range:
load('MatsRandomRes_4Neurons_4Paper_Rescaled_Sims.mat','results_rescaled');
results1 = results_rescaled;    clear results_rescaled
load('MatsRandomRes_4Neurons_4Paper_Rescaled_Sims_2.mat','results_rescaled');
results2 = results_rescaled';    clear results_rescaled
load('MatsRandomRes_4Neurons_4Paper.mat','results');
results3 = results;    clear results
results = [results1,results2,results3];  clear results1 results2 results3


c4p = class4paper(results);



