
clear all; close all; clc

% set default options
set(0,'defaultlinelinewidth',2);

% Load data:
GA1 = load('VGAM_08_01_08_00_GA_only.mat');
GA2 = load('VGAM_07_31_08_02_NN_only.mat');
GA3 = load('VGAM_07_30_20_59_rescale_only.mat');
GA4 = load('VGAM_07_30_08_43_NN_and_rescale.mat');

legends = {'MOGA','MOGA + NN',...
        'MOGA + re-scaling','MOGA + NN + re-scaling'};

    % the order of the fitnesses from the GA run:
fitnessOrder = {'VelFit','NrgEffFit',...
    'VelRangeFit #1','VelRangeFit #2','VelRangeFit #3',...
    'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
    'VelRangeFit #7','VelRangeFit #8','EigenFit'};
    
% which fit to plot:
FitNum = 3;
% get x-axis data:
x_data = 1:GA1.GA.Generations;

for i=1:11
    FitNum = i;
    % get y-axis data:
    y_data1 = squeeze(max(GA1.GA.Fit(:,FitNum,:),[],1));
    y_data2 = squeeze(max(GA2.GA.Fit(:,FitNum,:),[],1));
    y_data3 = squeeze(max(GA3.GA.Fit(:,FitNum,:),[],1));
    y_data4 = squeeze(max(GA4.GA.Fit(:,FitNum,:),[],1));

    % y_data1_mean = squeeze(mean(GA1.GA.Fit(:,3,:),1));
    % y_data2_mean = squeeze(mean(GA2.GA.Fit(:,3,:),1));
    % y_data3_mean = squeeze(mean(GA3.GA.Fit(:,3,:),1));
    % y_data4_mean = squeeze(mean(GA4.GA.Fit(:,3,:),1));

    % get y-data mean:
    y_data1_mean = mean(y_data1,2);
    y_data2_mean = mean(y_data2,2);
    y_data3_mean = mean(y_data3,2);
    y_data4_mean = mean(y_data4,2);

    % get y-data stDev:
    y_data1_std = std(y_data1,[],2);
    y_data2_std = std(y_data2,[],2);
    y_data3_std = std(y_data3,[],2);
    y_data4_std = std(y_data4,[],2);

    figure
    hold on
    h1=plot(x_data, y_data1);
    h2=plot(x_data, y_data2);
    h3=plot(x_data, y_data3);
    h4=plot(x_data, y_data4);
    grid minor;
    legend([h1,h2,h3,h4], legends, 'Location', 'Southeast');
    xlabel('Generation');   ylabel('Vel fitness');
    title(['Fitness "',fitnessOrder{1,i},'" over Generation']);
    set(gca,'FontSize',12, 'FontWeight','bold');
end


%% MOOGA figures, show divercity in clustering:
% https://www.mathworks.com/help/stats/kmeans.html

seq_1st_gen = []; 
seq_last_gen = GA4.GA.Seqs(:,:,end);

num_of_clusters = 4;
X = seq_last_gen;

param2plot = [1,2]; % according to seqOrder
Param1 = X(:,param2plot(1));
Param2 = X(:,param2plot(2));
XParam = '\tau';
YParam = 'b';

opts = statset('Display','final');
[idx,C,sumd,D] = kmeans([Param1,Param2],num_of_clusters,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);
% [idx,C,sumd,D] = kmeans(X,num_of_clusters,'Distance','sqeuclidean',...
%     'Replicates',5,'Options',opts);

figure;
plot(Param1,Param2,'k*','MarkerSize',5);
title(['clustering with ',num2str(num_of_clusters),' clusters']);
xlabel(XParam);
ylabel(YParam);

%% plot clusters centers:
figure; hold on;
colors = {'r.','b.','g.','m.','c.','k.'};
colors1 = {'r','b','g','m','c','k'};
legends = cell(1,num_of_clusters*2);

for i=1:num_of_clusters
    plot(Param1(idx==i,1),Param2(idx==i,1),...
        colors{1,i},'MarkerSize',12)
    legends{1,2*i-1} = sprintf('Cluster %d',i);
    legends{1,2*i} = sprintf('Centroids %d',i);
    
    plot(C(i,1),C(i,2),'ko',...
     'MarkerSize',10,'LineWidth',4,...
     'MarkerFaceColor',colors1{1,i},'MarkerEdgeColor','k');
end
grid minor;
title(['clustering with ',num2str(num_of_clusters),' clusters']);
xlabel(XParam);
ylabel(YParam);
legend(legends{1,:});

%% Color the clusters areas:

%Use kmeans to compute the distance from each centroid to
% points on a grid. To do this, pass the centroids (C) and points on
% a grid to kmeans, and implement one iteration of the algorithm.

p1 = min(Param1):0.001:max(Param1);
p2 = min(Param2):0.01:max(Param2);
[p1G,p2G] = meshgrid(p1,p2);
P_Grid = [p1G(:),p2G(:)]; % Defines a fine grid on the plot

idx2Region = kmeans(P_Grid,num_of_clusters,...
    'MaxIter',1,'Start',C);
    % Assigns each node in the grid to the closest centroid
    
% kmeans displays a warning stating that the algorithm did not
% converge, which you should expect since the
% software only implemented one iteration.    

figure; hold on;
for i=1:num_of_clusters
    plot(P_Grid(idx2Region==i,1),P_Grid(idx2Region==i,2),colors{1,i});
end
plot(Param1,Param2,'k*','MarkerSize',5);
xlabel(XParam);
ylabel(YParam);
legend(legends{1,1:2:end});
hold off;