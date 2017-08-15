clear all; close all; clc

colors = {'r.','b.','g.','m.','c.','k.'};
colors1 = {'r','b','g','m','c','k'};

% set default options
set(0,'defaultlinelinewidth',2);

% Load data:
% GA1 = load('VGAM_08_01_08_00_GA_only.mat');
% GA2 = load('VGAM_07_31_08_02_NN_only.mat');
% GA3 = load('VGAM_07_30_20_59_rescale_only.mat');
% GA4 = load('VGAM_07_30_08_43_NN_and_rescale.mat');

GA1 = load('VGAM_08_07_18_43_GA_only_NEW.mat');
GA2 = load('VGAM_08_11_10_50_NN_only_NEW.mat');
GA3 = load('VGAM_08_13_12_40_rescale_only_NEW.mat');
GA4 = load('VGAM_08_09_19_43_NN_and_rescale_NEW.mat');

legends = {'MOGA','MOGA + NN',...
        'MOGA + re-scaling','MOGA + NN + re-scaling'};
titleAdd = {'GA only','GA + NN','GA + rescale','everything'};

    % the order of the fitnesses from the GA run:
fitnessOrder = {'VelFit','NrgEffFit',...
    'VelRangeFit #1','VelRangeFit #2','VelRangeFit #3',...
    'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
    'VelRangeFit #7','VelRangeFit #8','EigenFit'};

seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};

seqOrder_extend = {'tau','b','c_1','c_2','c_3','c_4',...
                 'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
                 'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}',...
                 'ks_\tau','ks_c1','ks_c2','ks_c3','ks_c4'};

 % define cell array with all of the data:
chosen_method = {GA1,GA2,GA3,GA4};
             
% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;

% which fit to plot:
FitNum = 3;
% get x-axis data:
last_gen = 20;
x_data = 1:last_gen;

clear GA1 GA2 GA3 GA4
%% plot max fit over generation:
whichFit2Plot = 1:11;
plot_max_fit_over_gen(whichFit2Plot,fitnessOrder,...
    legends,...
    {GA1,GA2,GA3,GA4},...
    last_gen);

%% plot max and Mean fir over generation num:
whichFit2Plot = 1:3;
plot_mean_fit_over_gen(whichFit2Plot,fitnessOrder,...
    {GA1,GA2,GA3,GA4},...
    last_gen);

%% Plot specific Method parameters:
whichMethod = GA1;
titleCase = titleAdd{1,1};
normParam = param_norm_minMax(MML,whichMethod ,last_gen);

% Choose two parameters to plot:
Param_names = seqOrder_extend;
Param_name2Plot = {'tau','b'};
% get the parameters names:
X1name = Param_name2Plot{1,1};
X2name = Param_name2Plot{1,2};
% normalize the parameters between '0' to '1':
X1 = normParam(:,strcmp(X1name,Param_names));
X2 = normParam(:,strcmp(X2name,Param_names));

figure;
plot(X1,X2,'k*','MarkerSize',5);
title('param distribution');
xlabel(X1name);
ylabel(X2name);
axis([0,1,0,1]);
grid minor;

% MOOGA clustering:
% https://www.mathworks.com/help/stats/kmeans.html

num_of_clusters = 5;

opts = statset('Display','final');
[idx,C,sumd,D] = kmeans([X1,X2],num_of_clusters,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);

Tend_ratio = whichMethod.GA.Tend_ratio(:,1,last_gen);
% define "good" CPGs as one that T(end)/Sim.Tend is bigger than 0.99
ids_good = (Tend_ratio >0.99);

figure; hold on;
% plot clusters centers:
kMean_plot_clusters_and_centers(gca,X1,X2,...
    X1name,X2name,C,idx,ids_good,titleCase);


%% display parameters divercity:
num_of_clusters = 5;

figure; hold on;

Param_name2Plot = {'tau','b'};
X1name = Param_name2Plot{1,1};
X2name = Param_name2Plot{1,2};

for j=1:4
    
    allSeqs = chosen_method{1,j}.GA.Seqs;
    Param = zeros(size(allSeqs,1),size(allSeqs,2));
    
    % normalize the seq parameters to be between '0' to '1':
    normParam = param_norm_minMax(MML,chosen_method{1,j},last_gen);
    
    % extract the sequences of the wanted generation:
    Xtemp = chosen_method{1,j}.GA.Seqs(:,:,last_gen);
    
    % get the norm seq:
    X1 = normParam(:,strcmp(X1name,seqOrder_extend));
    X2 = normParam(:,strcmp(X2name,seqOrder_extend));

    % clustering using 'kmeans':
    opts = statset('Display','final');
    [idx,C,~,~] = kmeans([X1,X2],num_of_clusters,'Distance','sqeuclidean',...
        'Replicates',5,'Options',opts);
    
    % get "good CPGs" (CBs that walking):
    Tend_ratio = chosen_method{1,j}.GA.Tend_ratio(:,1,last_gen);
    ids_good = (Tend_ratio >0.99);
    
    ax = subplot(2,2,j); hold on;
    X_label = [X1name,' norm'];
    Y_label = [X2name,' norm'];
    kMean_plot_clusters_and_centers(ax,X1,X2,...
        X_label,Y_label,C,idx,ids_good,titleAdd{1,j});
    axis([0,1,0,1]);
    

end

%% display fitness divercity:
num_of_clusters = 5;

figure; hold on;

Param_name2Plot = {'VelFit','VelRangeFit #1'};

for j=1:4
    
    % exstract the fitnesses:
    X = chosen_method{1,j}.GA.Fit(:,:,last_gen);
    
    X1 = X(:,strcmp(Param_name2Plot{1,1},fitnessOrder));
    X2 = X(:,strcmp(Param_name2Plot{1,2},fitnessOrder));
    X1name = Param_name2Plot{1,1};
    X2name = Param_name2Plot{1,2};
    
    % clustering using 'kmeans':
    opts = statset('Display','final');
    [idx,C,~,~] = kmeans([X1,X2],num_of_clusters,'Distance','sqeuclidean',...
        'Replicates',5,'Options',opts);

    % get "good CPGs" (CBs that walking):
    Tend_ratio = chosen_method{1,j}.GA.Tend_ratio(:,1,last_gen);
    ids_good = (Tend_ratio >0.99);
    
    ax = subplot(2,2,j); hold on;
    X_label = X1name;
    Y_label = X2name;
    kMean_plot_clusters_and_centers(ax,X1,X2,...
        X_label,Y_label,C,idx,ids_good,titleAdd{1,j})
    axis(ax,[0,1,0,1]);

end

%% plot Pareto:
fit1Num = 1; % 'VelFit'
fit2Num = 3; % 'VelRangeFit#1'
MOGA_pareto_plot(chosen_method,fit1Num,fit2Num,...
    fitnessOrder,last_gen,titleAdd);

%% plot Tend ratio:
plot_Tend_ratio(chosen_method,last_gen,titleAdd);

%% plot generation runTime:
plot_gen_runTime(chosen_method,last_gen,titleAdd);

%% plot Hist of parameter over generation
whichParam = 'VelRangeFit #1';
plot_param_hist_over_genNum(chosen_method,whichParam,last_gen,...
    fitnessOrder,seqOrder_extend,titleAdd);
