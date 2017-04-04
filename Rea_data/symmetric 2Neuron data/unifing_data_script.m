
clear all;clc
load('MatsRandomRes_2Neurons_symm_trainData_wide_range1.mat')
results1=results; clear results
load('MatsRandomRes_2Neurons_symm_trainData_wide_range_all.mat')
results2=results; clear results

results = horzcat(results1,results2);
clear results1 results2

save('MatsRandomRes_2Neurons_symm_trainData_wide_range_all.mat','results')