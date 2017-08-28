clear all; close all; clc

% set default options
set(0,'defaultlinelinewidth',2);

InFiles_names = {'VGAM_08_17_17_41_GA_only_2000genes_10Gen.mat',...
    'VGAM_08_20_11_12_NN_only_2000genes_10Gen.mat',...
    'VGAM_08_17_17_44_rescale_only_2000genes_10Gen.mat',...
    'VGAM_08_22_22_36_NN_and_rescale_2000genes_10Gen.mat'};

GA_graphs = plotMOOGA4Paper(InFiles_names);

% which fit to plot:
FitNum = 3;
% get x-axis data:
last_gen = 10;
x_data = 1:last_gen;

% home many clusters in divesity plots:
num_of_clusters = 4;

%% plot max fit over generation:
whichFit2Plot = 1:11;
GA_graphs.plot_fit_over_gen(whichFit2Plot,last_gen)

%% plot max and Mean fit over generation num:
whichFit2Plot = 1:3;
GA_graphs.plot_mean_fit_over_gen(whichFit2Plot,last_gen,'all')

GA_graphs.plot_mean_fit_over_gen(whichFit2Plot,last_gen,'top_pop')
%% display parameters diversity:
whichParam = {'tau','b'};
GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,'param','plot by Tend ratio')
%% display parameters diversity for top 15%: 
whichParam = {'tau','b'};
GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,'param','plot by TopPop')

whichParam = {'VelFit','VelRangeFit #1'};
GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,'fit','plot by TopPop')

GA_graphs.map_divercity_based_on_topPop(last_gen,num_of_clusters,...
    {'tau','b'},'param',{'VelFit','VelRangeFit #1'},'fit')

GA_graphs.map_divercity_based_on_topPop(last_gen,num_of_clusters,...
    {'VelFit','VelRangeFit #1'},'fit',{'tau','b'},'param')

%% diversity 3D plot:
whichParam = {'VelFit','VelRangeFit #1','NrgEffFit'};
GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,'fit','plot by TopPop')

GA_graphs.map_divercity_based_on_topPop(last_gen,num_of_clusters,...
    whichParam,'fit',{'tau','b'},'param')

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
    whichParam,'param','plot by end condition')

whichParam = {'VelFit','VelRangeFit #1'};
GA_graphs.plot_divercity(last_gen,num_of_clusters,...
    whichParam,'fit','plot by end condition')

%% plot Pareto:
fit1Num = 1; % 'VelFit'
fit2Num = 3; % 'VelRangeFit#1'
GA_graphs.MOGA_pareto_plot(fit1Num,fit2Num,last_gen,'showAll');

GA_graphs.MOGA_pareto_plot(fit1Num,fit2Num,last_gen,'showBest');
%% plot Tend ratio:
GA_graphs.plot_Tend_ratio(last_gen);

%% plot generation runTime:
GA_graphs.plot_gen_runTime(last_gen);

%% plot Hist of parameter over generation
whichParam = 'VelRangeFit #1';
GA_graphs.plot_param_hist_over_genNum(whichParam,last_gen);


