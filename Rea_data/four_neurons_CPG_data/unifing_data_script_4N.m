
clear all;clc
load('MatsRandomRes_02_02_2017.mat')
results1=results; clear results
load('MatsRandomRes_03_02_2017.mat')
results2=results; clear results
load('MatsRandomRes_04_02_2017.mat')
results3=results; clear results

results = horzcat(results1,results2,results3);
clear results1 results2 results3

filename = 'MatsRandomRes_2-4_02_2017.mat'
save(filename,'results')