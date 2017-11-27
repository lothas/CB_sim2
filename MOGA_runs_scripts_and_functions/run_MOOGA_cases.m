
close all; clear all; clc

%%
% fileIn = 'VGAM_2N_symm_10_01_09_50_base_population_5000genes.mat';


% 
% % GA_try_2N_Symm_Matsuoka('GA + rescale',[]);
% GA_try_2N_Symm_Matsuoka('GA + NN_classi',[]);
% % GA_try_2N_Symm_Matsuoka('GA + NN_classi + rescale',[]);

%% 2N general CPG
% generate_GenomeFile('2N_general');

% GA_try_2N_General_Matsuoka('GA only',[]);
% GA_try_2N_General_Matsuoka('GA + NN_classi',[]);

%% 4N TagaLike CPG
% generate_GenomeFile('4N_tagaLike_1symmInput');
% 
% GA_try_TagaLike_Matsuoka('GA only',[]);
% GA_try_TagaLike_Matsuoka('GA only',[]);
% GA_try_TagaLike_Matsuoka('GA + NN_classi',[]);
% GA_try_TagaLike_Matsuoka('GA + NN_classi',[]);

%% 4N TagaLike CPG
generate_GenomeFile('4N_tagaLike_generalInput');
GA_try_TagaLike_Matsuoka('GA only',[]);
GA_try_TagaLike_Matsuoka('GA + NN_classi',[]);
% GA_try_TagaLike_Matsuoka('GA only',[]);
% GA_try_TagaLike_Matsuoka('GA + NN_classi',[]);

%% 6N TagaLike CPG
clc
generate_GenomeFile('6N_tagaLike_2Ank_torques');

trainDataFile = 'MatsRandomRes_6N_TagaLike_TrainingSet.mat';
GA_try_6N_TagaLike_Matsuoka('GA only',[],trainDataFile);
GA_try_6N_TagaLike_Matsuoka('GA + NN_classi',[],trainDataFile);
GA_try_6N_TagaLike_Matsuoka('GA only',[],trainDataFile);
GA_try_6N_TagaLike_Matsuoka('GA + NN_classi',[],trainDataFile);

%% 6N TagaLike CPG_ with semi-symm ankle/hip coupling (ensure phase diff) 
clc
generate_GenomeFile('6N_tagaLike_2Ank_torques_symm');
% trainDataFile = 'MatsRandomRes_6N_TagaLike_TrainingSet_test';
trainDataFile = 'MatsRandomRes_6N_TagaLike_TrainingSet_2';
GA_try_6N_TagaLike_Matsuoka('GA only',[],trainDataFile);
GA_try_6N_TagaLike_Matsuoka('GA + NN_classi',[],trainDataFile);
GA_try_6N_TagaLike_Matsuoka('GA only',[],trainDataFile);
GA_try_6N_TagaLike_Matsuoka('GA + NN_classi',[],trainDataFile);
% GA_try_6N_TagaLike_Matsuoka('GA only',[],trainDataFile);
% GA_try_6N_TagaLike_Matsuoka('GA + NN_classi',[],trainDataFile);
% GA_try_6N_TagaLike_Matsuoka('GA only',[],trainDataFile);
% GA_try_6N_TagaLike_Matsuoka('GA + NN_classi',[],trainDataFile);
% GA_try_6N_TagaLike_Matsuoka('GA only',[],trainDataFile);
% GA_try_6N_TagaLike_Matsuoka('GA + NN_classi',[],trainDataFile);