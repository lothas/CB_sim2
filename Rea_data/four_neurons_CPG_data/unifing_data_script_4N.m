
clear all;clc
load('MatsRandomRes_30_01_2017.mat')
results1=results; clear results
load('MatsRandomRes_31_01_2017.mat')
results2=results; clear results
load('MatsRandomRes_02_01_2017.mat')
results3=results; clear results

results = horzcat(results1,results2,results3);
clear results1 results2 results3

filename = 'MatsRandomRes_1-31_01_2017.mat';
save(filename,'results');