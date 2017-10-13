
close all; clear all; clc

%%
% fileIn = 'VGAM_2N_symm_10_01_09_50_base_population_5000genes.mat';



% GA_try_2N_Symm_Matsuoka('GA + rescale',[]);
GA_try_2N_Symm_Matsuoka('GA + NN_classi',[]);
% GA_try_2N_Symm_Matsuoka('GA + NN_classi + rescale',[]);

%% 2N general CPG
generate_GenomeFile('2N_general');

GA_try_2N_General_Matsuoka('GA only',[]);
GA_try_2N_General_Matsuoka('GA + NN_classi',[]);

