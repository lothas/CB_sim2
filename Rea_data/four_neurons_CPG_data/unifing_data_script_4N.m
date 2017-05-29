
clear all;clc
load('MatsRandomRes_1-31_01_2017.mat')
results1=results; clear results
load('MatsRandomRes_2-4_02_2017.mat')
results2=results; clear results
% load('MatsRandomRes_02_01_2017.mat')
% results3=results; clear results

results = horzcat(results1,results2);
clear results1 results2

filename = 'MatsRandomRes_all_from_1-2_2017.mat';
save(filename,'results');

%%

results2 = myCode.sim_results(myCode.ids);