
% % % Generating the genome:
% define Mutation strength:
MutDelta0 = 0.04;   MutDelta1 = 0.02;

nAnkle = 1;%1; % Number of ankle torques
nHip = 1;   % Number of hip torques
N = nAnkle+nHip;

genome_file = 'MatsuokaGenome_4Neuron_tagaLike.mat';
maxAnkle = 20;   % Max ankle torque
maxHip = 8;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;

mw = 0*ones(1,4);
Mw = 10*ones(1,4);

% CPG strucute: (ALSO Symm W_ij = W_ji)
%   H_F   H_E           % 
% 4 O-----O 3           %   
%    \    |             %   w = [0  , W12, 0  , 0  ; 
%     \   |             %        W21, 0  , w23, W24;
%      \  |             %        0  , 0  , 0  , W34;
%       \ |             %        0  , 0  , w43, 0  ;
% 1 O-----O 2           % w12=w21 = w1  
%  A_F    A_E           % w23 = w2
%                       % w24 = w3
%                       % w43=w34 = w4

% Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
Keys = {'\tau_r', 'beta',      'amp_same4each_joint',   '4neuron_taga_like', 'ks_\tau','ks_c_same4each_joint', 'IC_matsuoka';
              1 ,      1,                          2,                     4,        1 ,                     2,            0 };
Range = {  0.02 ,    0.2,        0*[maxAnkle,maxHip],                    mw,      -10 , -0.1*[maxAnkle,maxHip]; % Min
           0.25 ,    2.5,          [maxAnkle,maxHip],                    Mw,       10 ,  0.1*[maxAnkle,maxHip]}; % Max

       
save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');   

clear all
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

% % % LOAD files:
InFiles_names = {'VGAM_4N_TagaLike_10_29_20_38_general_tonicInputs__GA_only.mat',...
    'VGAM_4N_TagaLike_10_29_23_55_general_tonicInputs__GA_only.mat',...
    'VGAM_4N_TagaLike_10_30_08_41_general_tonicInputs__NN_classi_only.mat',...
    'VGAM_4N_TagaLike_10_30_13_43_general_tonicInputs__NN_classi_only.mat'};
Legends = {'GA1','GA2','NN1','NN2'};

% % % the order of the parametrs in CPG Sequence:
seqOrder = {'tau' ,'b', 'c1','c2', 'w1', 'w2', 'w3', 'w4',...
    'k_tau','k_{c1}','k_{c2}'};

GA_graphs = plotMOOGA4Paper(MML,InFiles_names,Legends,seqOrder);

% which fit to plot:
FitNum = 3;% get x-axis data:
last_gen = 10;
x_data = 1:last_gen;

% home many clusters in divesity plots:
num_of_clusters = 4;

%% 
GA_file_num = 2;
genNum = 10;
% GA_graphs.plot_seqs_in_gen(GA_file_num,genNum,1:16)

duration = 20;
geneID = 2;
% GA_graphs.animate_seq(GA_file_num,genNum, geneID, duration, [])

clear GA_file_num genNum duration geneID
%% plot max fit over generation:
close all
whichFit2Plot = 3;%1:11;
GA_graphs.plot_fit_over_gen(whichFit2Plot,last_gen);

%% plot max and Mean fit over generation num:
close all
whichFit2Plot = 3;%1:3;
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

% % % the order of the parametrs in CPG Sequence:
% Symm tonic inputs:
seqOrder = {'tau' ,'b', 'c1','c2', 'w1', 'w2', 'w3', 'w4',...
    'k_tau','k_{c1}','k_{c2}'};

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;

% % LOAD GA_only files:
% InFiles_names = {'VGAM_4N_TagaLike_10_24_10_25_1tonicInput__GA_only.mat',...
%     'VGAM_4N_TagaLike_10_23_17_27_1tonicInput__GA_only.mat',...
%     'VGAM_4N_TagaLike_10_24_00_32_1tonicInput__GA_only.mat'};
% Legends = {'GA1','GA2','GA3'};
% LOAD GA_only files:
InFiles_names = {'VGAM_4N_TagaLike_10_24_10_25_1tonicInput__GA_only.mat',...
    'VGAM_4N_TagaLike_10_23_17_27_1tonicInput__GA_only.mat'};
Legends = {'GA1','GA2'};
GA_only = plotMOOGA4Paper(MML,InFiles_names,Legends,seqOrder);

% % LOAD GA_NN files:
InFiles_names = {'VGAM_4N_TagaLike_10_24_10_15_1tonicInput__NN_classi_only_wrong_C.mat',...
    'VGAM_4N_TagaLike_10_24_14_29_1tonicInput__NN_classi_only_wrong_C.mat',...
    'VGAM_4N_TagaLike_10_24_20_35_1tonicInput__NN_classi_only_wrong_C.mat'};
Legends = {'NN1','NN2','NN3'};
GA_NN = plotMOOGA4Paper(MML,InFiles_names,Legends,seqOrder);

% get x-axis data:
last_gen = 10;

close all;
Title = {'velFit','NrgFit','rangeVelFit'};

% % Plot mean of maximum fitness for all MOGA runs:
for j=1:3
    whichFit2Plot = j;
    [x_data,GA_only_max] = GA_only.plot_fit_over_gen(whichFit2Plot,last_gen);
    [~,GA_NN_max] = GA_NN.plot_fit_over_gen(whichFit2Plot,last_gen);

    GA_only_mean = mean(GA_only_max);
    GA_only_std = std(GA_only_max);

    GA_NN_mean = mean(GA_NN_max);
    GA_NN_std = std(GA_NN_max);

    figure;
    plot(x_data,GA_only_mean); hold on;
%     errorbar(x_data,GA_only_mean,GA_only_std);
    plot(x_data,GA_NN_mean);
%     errorbar(x_data,GA_NN_mean,GA_NN_std);
    xlabel('Generation number');
    ylabel('Velocity Fitness Max');
    legend('MOGA','MOGA with NN assist');
%     title([Title{1,j},'_{mean}']);
    grid minor
    axis([1,10,0,0.4]);
    set(gca,'fontsize',13)
    
    clear x_data GA_only_max GA_NN_max 
    clear GA_only_mean GA_only_std GA_NN_mean GA_NN_std
end


% % Plot mean of "TopPop Mean Fitness" (TPMF) fitness for all MOGA runs:
if false   % Don't plot it right now
    for j=1:3
        whichFit2Plot = j;
        [x_data,GA_only_TPMF] = ...
            GA_only.plot_mean_fit_over_gen(whichFit2Plot,last_gen,'top_pop');
        [~,GA_NN_TPMF] = ...
            GA_NN.plot_mean_fit_over_gen(whichFit2Plot,last_gen,'top_pop');

        GA_only_mean = mean(GA_only_TPMF);
        GA_only_std = std(GA_only_TPMF);

        GA_NN_mean = mean(GA_NN_TPMF);
        GA_NN_std = std(GA_NN_TPMF);

        figure;
        plot(x_data,GA_only_mean); hold on;
        plot(x_data,GA_NN_mean);
        xlabel('gen num');
        ylabel('mean of TopPop mean fitness');
        legend('GA only','NN');
        title([Title{1,j},'_{mean}']);
        grid minor

        clear x_data GA_only_TPMF GA_NN_TPMF 
        clear GA_only_mean GA_only_std GA_NN_mean GA_NN_std
    end
end