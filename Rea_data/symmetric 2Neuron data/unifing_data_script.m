
clear all;clc
load('MatsRandomRes_2Neurons_symm_test_for_FFT_new_1.mat')
results1=results; clear results
load('MatsRandomRes_2Neurons_symm_test_for_FFT_new_2.mat')
results2=results; clear results
load('MatsRandomRes_2Neurons_symm_test_for_FFT_new_3.mat')
results3=results; clear results

results = horzcat(results1,results2,results3);
clear results1 results2 results3

save('MatsRandomRes_2Neurons_symm_test_for_FFT_new_all.mat','results')