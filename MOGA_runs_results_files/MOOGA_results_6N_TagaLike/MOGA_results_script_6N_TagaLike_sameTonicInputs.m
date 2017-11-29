
% % % Generating the genome:
% define Mutation strength:
MutDelta0 = 0.04;   MutDelta1 = 0.02;

nAnkle = 2; % Number of ankle torques
nHip = 1;   % Number of hip torques
N = nAnkle+nHip;

genome_file = 'MatsuokaGenome_4Neuron_tagaLike.mat';
maxAnkle = 20;   % Max ankle torque
maxHip = 8;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;

mw = 0*ones(1,6);
Mw = 5*ones(1,6);

% CPG strucute: (ALSO Symm W_ij = W_ji)
% 1) every pair of Extensor and reflexor neurons are connected to
%   each other with symmetric weights.
% 2) both of the ankles' Extensor neurons are connected to both of the hip neurons
%
%  A2_F    A2_E
%(3) O-----O (4)
%       /  |
%      /   |
%     /    |
%    /     |
%   H_F   H_E
%(5) O-----O (6)  
%     \    |             %   w = [0  , W12, 0  , 0  , 0  , 0; 
%      \   |             %        W21, 0  , 0  , 0  , W25, W26;
%       \  |             %        0  , 0  , 0  , W34, 0  , 0;
%        \ |             %        0  , 0  , W43, 0  , W45, W46;
%(1) O-----O (2)         %        0  , 0  , 0  , 0  , 0  , W56;
%   A1_F    A1_E         %        0  , 0  , 0  , 0  , W65, 0;
%                        
%                        w12 = w21 = w34 = w43 = w1  
%                        w56 = 65  = w2
%                        w25 = w3
%                        w26 = w4
%                        w45 = w5
%                        w46 = w6

% Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
Keys = {'\tau_r', 'beta', 'amp_6n_symm',   '6neuron_taga_like', 'ks_\tau',     'ks_c_6n_symm', 'IC_matsuoka';
              1 ,      1,             1,                     6,        1 ,                 1 ,            0 };
Range = {  0.02 ,    0.2,             0,                    mw,      -10 ,                 -0.1*maxAnkle; % Min
           0.25 ,    2.5,      maxAnkle,                    Mw,       10 ,                  0.1*maxAnkle}; % Max

       
save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');   

clear all
%%
clear all; close all; clc

% set default options
set(0,'defaultlinelinewidth',2);

load('MatsuokaGenome_4Neuron_tagaLike.mat','N');

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 2*N;
clear N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % LOAD files:
% InFiles_names = {'VGAM_6N_TagaLike_11_09_17_01_same_tonicInputs__GA_only.mat',...
%     'VGAM_6N_TagaLike_11_10_00_50_same_tonicInputs__GA_only.mat',...
%     'VGAM_6N_TagaLike_11_10_15_51_same_tonicInputs__GA_only.mat',...
%     'VGAM_6N_TagaLike_11_11_00_41_same_tonicInputs__NN_classi_only.mat',...
%     'VGAM_6N_TagaLike_11_11_12_43_same_tonicInputs__NN_classi_only.mat',...
%     'VGAM_6N_TagaLike_11_11_21_58_same_tonicInputs__NN_classi_only.mat'};
% Legends = {'GA1','GA2','GA3','NN1','NN2','NN3'};
% 
% seqOrder = {'tau' ,'b', 'c', 'w1', 'w2', 'w3', 'w4', 'w5', 'w6','k_tau','k_{c}'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% InFiles_names = {'VGAM_6N_TagaLike_11_12_14_11_same_tonicInputs__GA_only.mat',...
%     'VGAM_6N_TagaLike_11_12_17_47_same_tonicInputs__NN_classi_only.mat'};
% Legends = {'GA1 test','NN1 test'};

% %
% LOAD GA_only files:
InFiles_names = {'VGAM_6N_TagaLike_11_15_08_15_same_tonicInputs_20Gen_500Genes__GA_only.mat',...
    'VGAM_6N_TagaLike_11_15_18_46_same_tonicInputs_20Gen_500Genes__GA_only.mat',...
    'VGAM_6N_TagaLike_11_16_04_44_same_tonicInputs_20Gen_500Genes__GA_only.mat',...
    'VGAM_6N_TagaLike_11_22_17_20_same_tonicInputs_20Gen_500Genes__GA_only.mat',...
    'VGAM_6N_TagaLike_11_23_02_19_same_tonicInputs_20Gen_500Genes__GA_only.mat',...
    'VGAM_6N_TagaLike_11_23_13_40_same_tonicInputs_20Gen_500Genes__GA_only.mat',...
    'VGAM_6N_TagaLike_11_23_23_48_same_tonicInputs_20Gen_500Genes__GA_only.mat',...
    'VGAM_6N_TagaLike_11_24_11_16_same_tonicInputs_20Gen_500Genes__GA_only.mat',...
    'VGAM_6N_TagaLike_11_25_19_42_same_tonicInputs_20Gen_500Genes__GA_only.mat'};
Legends = {'GA1','GA2','GA3','GA4','GA5','GA6','GA7','GA8','GA9'};

% % % the order of the parametrs in CPG Sequence:
seqOrder = {'tau' ,'b', 'c', 'w1', 'w2', 'w3', 'w4','k_tau','k_{c}'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GA_graphs = plotMOOGA4Paper(MML,InFiles_names,Legends,seqOrder);

% which fit to plot:
FitNum = 3;% get x-axis data:
last_gen = 20;
x_data = 1:last_gen;

% home many clusters in divesity plots:
num_of_clusters = 4;

%% 
% close all
GA_file_num = 2;
genNum = 5;
GA_graphs.plot_seqs_in_gen(GA_file_num,genNum,4)

clear GA_file_num genNum
%%
clc
GA_file_num = 2;
genNum = 5;
duration = 10;
timestep = 0.05;
geneID = 87;
% simRun = GA_graphs.animate_seq(GA_file_num,genNum, geneID,duration,timestep, []);
% simRun = GA_graphs.stickFigure_plot(GA_file_num,genNum, geneID,duration,timestep, []);
simRun = GA_graphs.CBstick_Figure_plot(GA_file_num,genNum, geneID,duration,timestep, []);

%%
n_ankle = 2;
n_hip = 1;
N = n_ankle + n_hip;
Torque = cell(1,N);
TorquesNames = {'\tau_{ank1}','\tau_{ank2}','\tau_{hip}'};

for i=1:N
    Torque{1,i} = max(simRun.Out.X(:,4+(2*i-1)),0)-...
        max(simRun.Out.X(:,4+(2*i)),0);
end
T_ankle = max(simRun.Out.X(:,5),0)-max(simRun.Out.X(:,6),0);
T_hip = max(simRun.Out.X(:,7),0)-max(simRun.Out.X(:,8),0);

figure;
subplot(2,1,1);
plot(simRun.Out.T,simRun.Out.Torques(:,1)); hold on;
plot(simRun.Out.T,simRun.Out.Torques(:,2));
legend ('\tau_{\theta_1}','\tau_{\theta_2}');
xlabel('time'); ylabel('Torque');
grid minor;

subplot(2,1,2); hold on;
for i=1:N
    plot(simRun.Out.T,Torque{1,i});
end
legend (TorquesNames);
xlabel('time'); ylabel('Torque');
grid minor;

% clear GA_file_num genNum duration geneID
%% plot max fit over generation:
close all; clc
whichFit2Plot = 1:3;%1:11;
GA_graphs.plot_fit_over_gen(whichFit2Plot,last_gen);

%% plot max and Mean fit over generation num:
close all
whichFit2Plot = 1;%1:3;
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
