
clear all; close all; clc

colors = {'r.','b.','g.','m.','c.','k.'};
colors1 = {'r','b','g','m','c','k'};

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

seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};
    
% which fit to plot:
FitNum = 3;
% get x-axis data:
x_data = 1:GA1.GA.Generations;

if false
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

end

%% MOOGA figures, show divercity in clustering:
% https://www.mathworks.com/help/stats/kmeans.html

chosen_method = GA1;
Param = Fit_last_gen;%seq_last_gen;
Param_names = fitnessOrder;
Param_name2Plot = {'VelFit','NrgEffFit'};

%         % the order of the fitnesses from the GA run:
%          fitnessOrder = {'VelFit','NrgEffFit',...
%              'VelRangeFit #1','VelRangeFit #2','VelRangeFit #3',...
%              'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
%              'VelRangeFit #7','VelRangeFit #8','EigenFit'};
%          
%          seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
%              'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
%              'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};

seq_last_gen = chosen_method.GA.Seqs(:,:,end);
Fit_last_gen = chosen_method.GA.Fit(:,:,end);


num_of_clusters = 4;
X = Param;

Param1 = X(:,strcmp(Param_name2Plot{1,1},Param_names));
Param2 = X(:,strcmp(Param_name2Plot{1,2},Param_names));
XParam = Param_name2Plot{1,1};
YParam = Param_name2Plot{1,2};

opts = statset('Display','final');
[idx,C,sumd,D] = kmeans([Param1,Param2],num_of_clusters,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);

% figure;
% plot(Param1,Param2,'k*','MarkerSize',5);
% title(['clustering with ',num2str(num_of_clusters),' clusters']);
% xlabel(XParam);
% ylabel(YParam);

%% plot clusters centers:
figure; hold on;
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


%%

num_of_clusters = 3;


figure; hold on;
legends = cell(1,num_of_clusters*2);
tiltleAdd = {'GA only','GA + NN','GA + rescale','everything'};

chosen_method = {GA1,GA2,GA3,GA4};

method = 'Fit';
Param_names = fitnessOrder;
Param_name2Plot = {'VelFit','NrgEffFit'};

% method = 'Seq';
% Param_names = seqOrder;
% Param_name2Plot = {'tau','b'};

for j=1:4
    
    switch method
        case 'Seq'
            Param = chosen_method{1,j}.GA.Seqs(:,:,end);
        case 'Fit'
            Param = chosen_method{1,j}.GA.Fit(:,:,end);
    end
    
    X = Param;
    
    Param1 = X(:,strcmp(Param_name2Plot{1,1},Param_names));
    Param2 = X(:,strcmp(Param_name2Plot{1,2},Param_names));
    XParam = Param_name2Plot{1,1};
    YParam = Param_name2Plot{1,2};

    opts = statset('Display','final');
    [idx,C,~,~] = kmeans([Param1,Param2],num_of_clusters,'Distance','sqeuclidean',...
        'Replicates',5,'Options',opts);

    ax = subplot(2,2,j); hold on;
    for i=1:num_of_clusters
        plot(ax,Param1(idx==i,1),Param2(idx==i,1),...
            colors{1,i},'MarkerSize',12)
        legends{1,2*i-1} = sprintf('Cluster %d',i);
        legends{1,2*i} = sprintf('Centroids %d',i);

        plot(ax,C(i,1),C(i,2),'ko',...
         'MarkerSize',10,'LineWidth',4,...
         'MarkerFaceColor',colors1{1,i},'MarkerEdgeColor','k');
    end
    
    grid minor;
    title(['clustering with ',num2str(num_of_clusters),...
        ' clusters: ',tiltleAdd{1,j}]);
    xlabel(XParam);
    ylabel(YParam);
    legend(legends{1,:},'Location','best');

end


%% Competetive layers:
num_of_clusters = 3;

chosen_method = GA1;
Param_names = fitnessOrder;
Param_name2Plot = {'VelFit','NrgEffFit'};

% X = chosen_method.GA.Seqs(:,:,end);
X = chosen_method.GA.Fit(:,:,end);

Param1 = X(:,strcmp(Param_name2Plot{1,1},Param_names));
Param2 = X(:,strcmp(Param_name2Plot{1,2},Param_names));
XParam = Param_name2Plot{1,1};
YParam = Param_name2Plot{1,2};

inputs = ([Param1,Param1])';

net = competlayer(num_of_clusters);
net = train(net,inputs);
view(net)
outputs = net(inputs);
idx = vec2ind(outputs);

figure; hold on;
legends = cell(1,num_of_clusters*2);

for i=1:num_of_clusters
    plot(Param1(idx==i,1),Param2(idx==i,1),...
        colors{1,i},'MarkerSize',12)
    legends{1,2*i-1} = sprintf('Cluster %d',i);
    legends{1,2*i} = sprintf('Centroids %d',i);
end
grid minor;
title(['clustering with ',num2str(num_of_clusters),' clusters']);
xlabel(XParam);
ylabel(YParam);
legend(legends{1,1:2:end});
    