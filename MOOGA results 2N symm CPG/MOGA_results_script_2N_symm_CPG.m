

%% create Genome file if necessary:
genome_file = 'MatsuokaGenome_2Neuron_Symm.mat';
nAnkle = 1;%1; % Number of ankle torques
nHip = 1;   % Number of hip torques
maxAnkle = 10;   % Max ankle torque
maxHip = 10;    % Max hip torque
Mamp = [10,10];
mamp = 0*Mamp;
N = nAnkle+nHip;
% %     % 2neuron symmetric specific range%%

       % % % Narrow b Narrow W Narrow tau
Mw = 5*1;
mw = 0*Mw;
Keys = {'\tau_r', 'beta','amp_2n_same_inputs',    '2neuron_symm_weights', 'ks_\tau',     'ks_c_2n_symm', 'IC_matsuoka';
              1 ,      1,                   2,                         1,        1 ,          1,            0 };
Range = {  0.02 ,    0.2,               [0,0],                         0,   -0.001 ,       -0.2; % Min
           0.10  ,   2.5,             [10,10],                         5,    0.001 ,       0.2}; % Max

       
MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

%%
clear all; close all; clc

% set default options
set(0,'defaultlinelinewidth',2);

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;

% % % LOAD baseSeq comparison analisys:
% InFiles_names = {'VGAM_2N_symm_10_01_10_01_GA_only from base pop.mat',...
%     'VGAM_2N_symm_10_01_10_53_NN_classi_from_base_pop.mat',...
%     'VGAM_2N_symm_10_01_10_22_rescale_from_base_pop.mat',...
%     'VGAM_2N_symm_10_01_11_16_NN_classi_rescale_from_base_pop.mat'};
% Legends = {'GA only','NN classi','rescale','rescale+NN'};

% % LOAD files:
% InFiles_names = {'VGAM_2N_symm_09_28_14_24_GA_only.mat',...
%     'VGAM_2N_symm_09_29_18_44_NN_classi_only.mat',...
%     'VGAM_2N_symm_09_29_18_46_rescale_only.mat',...
%     'VGAM_2N_symm_09_30_20_08_NNclassi_and_rescale.mat'};
% Legends = {'GA only','NN classi','rescale','rescale+NN'};


% % LOAD GA_only files:
% InFiles_names = {'VGAM_2N_symm_09_27_15_23_GA_only.mat',...
%     'VGAM_2N_symm_09_28_09_05_GA_only.mat',...
%     'VGAM_2N_symm_09_28_14_24_GA_only.mat'};
% Legends = {'VGAM 2N symm 09/27 15:23 GA only',...
%     'VGAM 2N symm 09/28 09:05 GA only',...
%     'VGAM 2N symm 09/28 14:24 GA only'};

% % % LOAD GA_rescale files:
% InFiles_names = {'VGAM_2N_symm_09_27_15_26_rescale_only.mat',...
%     'VGAM_2N_symm_09_28_09_06_rescale_only.mat',...
%     'VGAM_2N_symm_09_28_21_30_rescale_only.mat',...
%     'VGAM_2N_symm_09_29_18_46_rescale_only.mat',...
%     'VGAM_2N_symm_09_30_20_12_rescale_only.mat'};
% Legends = {'VGAM_2N_symm_09_27_15_26_rescale_only',...
%     'VGAM_2N_symm_09_28_09_06_rescale_only',...
%     'VGAM_2N_symm_09_28_21_30_rescale_only',...
%     'VGAM_2N_symm_09_29_18_46_rescale_only',...
%     'VGAM_2N_symm_09_30_20_12_rescale_only'};

% % LOAD GA_NN files:
InFiles_names = {'VGAM_2N_symm_09_28_21_24_NN_classi_only.mat',...
    'VGAM_2N_symm_09_28_21_27_NN_classi_only.mat',...
    'VGAM_2N_symm_09_29_18_44_NN_classi_only.mat',...
    'VGAM_2N_symm_09_29_18_49_NN_classi_only.mat'};
Legends = {'VGAM_2N_symm_09_28_21_24_NN_classi_only',...
    'VGAM_2N_symm_09_28_21_27_NN_classi_only',...
    'VGAM_2N_symm_09_29_18_44_NN_classi_only',...
    'VGAM_2N_symm_09_29_18_49_NN_classi_only'};

% % LOAD GA_NN+reacale files:
% InFiles_names = {'VGAM_2N_symm_09_30_10_08_NNclassi_and_rescale.mat',...
%     'VGAM_2N_symm_09_30_20_08_NNclassi_and_rescale.mat'};
% Legends = {'VGAM_2N_symm_09_30_10_08_NNclassi_and_rescale',...
%     'VGAM_2N_symm_09_30_20_08_NNclassi_and_rescale'};

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau' ,'b', 'c', 'NR', 'a',...
    'k_tau','k_{c}'};
% "NR" - not relevnt param 

GA_graphs = plotMOOGA4Paper(MML,InFiles_names,Legends,seqOrder);

% which fit to plot:
FitNum = 3;
% get x-axis data:
last_gen = 20;
x_data = 1:last_gen;

% home many clusters in divesity plots:
num_of_clusters = 4;

%% plot max fit over generation:
close all
whichFit2Plot = 1:3;%1:11;
[x_data,y_data] = GA_graphs.plot_fit_over_gen(whichFit2Plot,last_gen);

%% plot max and Mean fit over generation num:
close all
whichFit2Plot = 1:3;
% GA_graphs.plot_mean_fit_over_gen(whichFit2Plot,last_gen,'all')

[x_data,y_data_mean] = ...
    GA_graphs.plot_mean_fit_over_gen(whichFit2Plot,last_gen,'top_pop');
%% display parameters diversity:
whichParam = {'tau','b'};
normID = [1,1]; % norm or not
GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,normID,'plot by Tend ratio');

clear normID whichParam
%% display parameters diversity for top 15%: 
whichParam = {'tau','b'};   normID = [1,1]; % norm or not
% whichParam = {'VelFit','VelRangeFit #1'};   normID = [0,0];
% whichParam = {'tau','VelRangeFit #1'};  normID = [1,0];

GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,normID,'plot by TopPop')

% GA_graphs.map_divercity_based_on_topPop(last_gen,num_of_clusters,...
%     {'tau','b'},[1,1],{'VelFit','VelRangeFit #1'},[0,0]);
% 
% GA_graphs.map_divercity_based_on_topPop(last_gen,num_of_clusters,...
%     {'VelFit','VelRangeFit #1'},[0,0],{'tau','b'},[1,1]);

clear normID whichParam
%% diversity 3D plot:
whichParam = {'VelFit','VelRangeFit #1','NrgEffFit'};
GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,[0,0,0],'plot by TopPop');

GA_graphs.map_divercity_based_on_topPop(last_gen,num_of_clusters,...
    whichParam,[0,0,0],{'tau','b'},[1,1]);

%% display divercity plot from endCond:
% every 'Sim' as a 'endCond' as follow:
% Out.Type:
            %   -1 - 'Simulation stopped by user'
            %   0 - 'Reached end of tspan'
            %   1 - 'Robot fell down (hip height too low)'
            %   2 - 'Step length too small'
            %   3 - 'Robot hit the ground before extending the leg'
            %   4 - 'Finished taking _ steps on slope of _ rad' (when EndCond = 1)
            %   5 - 'Reached steady state limit cycle of period _ after _ steps'
            %   6 - 'GO: Sim was converging' (check if necessary)
            %   7 - 'NO-GO: Simulation was diverging' (check if necessary)
            %   8 - 'ZMP crossed the limits'
            
% We are looking for (endCond == 0)
whichParam = {'tau','b'};
GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,[1,1],'plot by end condition');

whichParam = {'VelFit','VelRangeFit #1'};
GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,[0,0],'plot by end condition');

%% plot Pareto:
close all
fit1Num = 1; % 'VelFit'
fit2Num = 3; % 'VelRangeFit#1'
GA_graphs.MOGA_pareto_plot(fit1Num,fit2Num,last_gen,'showAll');

GA_graphs.MOGA_pareto_plot(fit1Num,fit2Num,last_gen,'showBest_subPlots');

GA_graphs.MOGA_pareto_plot(fit1Num,fit2Num,last_gen,'showBest_onOnePlot');

clear fit1Num fit2Num
%% plot Tend ratio:
close all
GA_graphs.plot_Tend_ratio(last_gen);

%% plot generation runTime:
GA_graphs.plot_gen_runTime(last_gen);

%% plot Hist of parameter over generation
close all

whichParam = 'b';
GA_graphs.plot_param_hist_over_genNum(whichParam,last_gen);

%% plot Hist of parameter over generation (ALL)
for i=1:numel(seqOrder)
    GA_graphs.plot_param_hist_over_genNum(seqOrder{1,i},last_gen);
end
clear i


%% Analize the statistical perf of many:
clear all; close all; clc

% set default options
set(0,'defaultlinelinewidth',2);

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau' ,'b', 'c', 'NR', 'a',...
    'k_tau','k_{c}'};
% "NR" - not relevnt param 


% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;

% LOAD GA_only files:
InFiles_names = {'VGAM_2N_symm_09_27_15_23_GA_only.mat',...
    'VGAM_2N_symm_09_28_09_05_GA_only.mat',...
    'VGAM_2N_symm_09_28_14_24_GA_only.mat'};
Legends = {'VGAM 2N symm 09/27 15:23 GA only',...
    'VGAM 2N symm 09/28 09:05 GA only',...
    'VGAM 2N symm 09/28 14:24 GA only'};
GA_only = plotMOOGA4Paper(MML,InFiles_names,Legends,seqOrder);

% % LOAD GA_rescale files:
InFiles_names = {'VGAM_2N_symm_09_27_15_26_rescale_only.mat',...
    'VGAM_2N_symm_09_28_09_06_rescale_only.mat',...
    'VGAM_2N_symm_09_28_21_30_rescale_only.mat',...
    'VGAM_2N_symm_09_29_18_46_rescale_only.mat',...
    'VGAM_2N_symm_09_30_20_12_rescale_only.mat'};
Legends = {'VGAM_2N_symm_09_27_15_26_rescale_only',...
    'VGAM 2N symm 09/28 09:06 rescale only',...
    'VGAM 2N symm 09/28 21:30 rescale only',...
    'VGAM 2N symm 09/29 18:46 rescale only',...
    'VGAM 2N symm 09/30 20:12 rescale only'};
GA_rescale = plotMOOGA4Paper(MML,InFiles_names,Legends,seqOrder);

% % LOAD GA_NN files:
InFiles_names = {'VGAM_2N_symm_09_28_21_24_NN_classi_only.mat',...
    'VGAM_2N_symm_09_28_21_27_NN_classi_only.mat',...
    'VGAM_2N_symm_09_29_18_44_NN_classi_only.mat',...
    'VGAM_2N_symm_09_29_18_49_NN_classi_only.mat'};
Legends = {'VGAM 2N symm 09/28 21:24 NN classi only',...
    'VGAM 2N symm 09/28 21:27 NN classi only',...
    'VGAM 2N symm 09/29 18:44 NN classi only',...
    'VGAM 2N symm 09/29 18:49 NN classi only'};
GA_NN = plotMOOGA4Paper(MML,InFiles_names,Legends,seqOrder);

% LOAD GA_NN+reacale files:
InFiles_names = {'VGAM_2N_symm_09_30_10_08_NNclassi_and_rescale.mat',...
    'VGAM_2N_symm_09_30_20_08_NNclassi_and_rescale.mat'};
Legends = {'VGAM 2N symm 09/30 10:08 NNclassi and rescale',...
    'VGAM 2N symm 09/30 20:08 NNclassi and rescale'};
GA_NN_rescale = plotMOOGA4Paper(MML,InFiles_names,Legends,seqOrder);

% get x-axis data:
last_gen = 20;

whichFit2Plot = 1:3;
[x_data,GA_only_max] = GA_only.plot_fit_over_gen(whichFit2Plot,last_gen);
[~,GA_rescale_max] = GA_rescale.plot_fit_over_gen(whichFit2Plot,last_gen);
[~,GA_NN_max] = GA_NN.plot_fit_over_gen(whichFit2Plot,last_gen);
[~,GA_NN_rescale_max] = GA_NN_rescale.plot_fit_over_gen(whichFit2Plot,last_gen);

close all;

GA_only_mean = mean(GA_only_max);
GA_only_std = std(GA_only_max);

GA_rescale_mean = mean(GA_rescale_max);
GA_rescale_std = std(GA_rescale_max);

GA_NN_mean = mean(GA_NN_max);
GA_NN_std = std(GA_NN_max);

GA_NN_rescale_mean = mean(GA_NN_rescale_max);
GA_NN_rescale_std = std(GA_NN_rescale_max);

figure;
plot(x_data,GA_only_mean); hold on;
plot(x_data,GA_rescale_mean);
plot(x_data,GA_NN_mean);
plot(x_data,GA_NN_rescale_mean);
xlabel('gen num');
ylabel('mean of Max fitness');
legend('GA only','rescale','NN','NN+rescale');