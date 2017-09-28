clear all; close all; clc

% set default options
set(0,'defaultlinelinewidth',2);

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;

InFiles_names = {'VGAM_2N_symm_09_27_15_23_GA_only.mat',...
    'VGAM_2N_symm_09_28_09_05_GA_only.mat',...
    'VGAM_2N_symm_09_27_15_28_NN_classi_only.mat',...
    'VGAM_2N_symm_09_27_15_26_rescale_only.mat'};

Legends = {'GA only 1','GA only 1','NN classi','rescale'};

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
GA_graphs.plot_fit_over_gen(whichFit2Plot,last_gen)

%% plot max and Mean fit over generation num:
close all
whichFit2Plot = 1:3;
% GA_graphs.plot_mean_fit_over_gen(whichFit2Plot,last_gen,'all')

GA_graphs.plot_mean_fit_over_gen(whichFit2Plot,last_gen,'top_pop')
%% display parameters diversity:
whichParam = {'tau','b'};
normID = [1,1]; % norm or not
GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,normID,'plot by Tend ratio');
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
