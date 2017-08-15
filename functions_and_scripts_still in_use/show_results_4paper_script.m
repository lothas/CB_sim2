
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

%%
close all;
if true
for i = 1:3 % get max fitness over generation:
    FitNum = i;
    % get y-axis data:
    y_data1 = squeeze(max(GA1.GA.Fit(:,FitNum,1:last_gen),[],1));
    y_data2 = squeeze(max(GA2.GA.Fit(:,FitNum,1:last_gen),[],1));
    y_data3 = squeeze(max(GA3.GA.Fit(:,FitNum,1:last_gen),[],1));
    y_data4 = squeeze(max(GA4.GA.Fit(:,FitNum,1:last_gen),[],1));

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

% for i=1:3 % get mean+std fitness over generation:
%     FitNum = i;
%     
%     fit_GA1 = GA1.GA.Fit(:,FitNum,1:last_gen);
%     fit_GA2 = GA2.GA.Fit(:,FitNum,1:last_gen);
%     fit_GA3 = GA3.GA.Fit(:,FitNum,1:last_gen);
%     fit_GA4 = GA4.GA.Fit(:,FitNum,1:last_gen);
%     
%     % get Max fit:
%     y_data1 = squeeze(max(fit_GA1,[],1));
%     y_data2 = squeeze(max(fit_GA2,[],1));
%     y_data3 = squeeze(max(fit_GA3,[],1));
%     y_data4 = squeeze(max(fit_GA4,[],1));
%     
%     % get Mean fit:
%     y_data1_mean = squeeze(mean(fit_GA1,1));
%     y_data2_mean = squeeze(mean(fit_GA2,1));
%     y_data3_mean = squeeze(mean(fit_GA3,1));
%     y_data4_mean = squeeze(mean(fit_GA4,1));
% 
%     % get STdev fit:
%     y_data1_std = squeeze(std(fit_GA1,[],1));
%     y_data2_std = squeeze(std(fit_GA2,[],1));
%     y_data3_std = squeeze(std(fit_GA3,[],1));
%     y_data4_std = squeeze(std(fit_GA4,[],1));
% 
%     figure; hold on;
% %     h1=plot(x_data, y_data1_mean,'b',x_data, y_data1_mean,'*b');
% %     h2=plot(x_data, y_data2_mean,'g',x_data, y_data2_mean,'*g');
% %     h3=plot(x_data, y_data3_mean,'r',x_data, y_data3_mean,'*r');
% %     h4=plot(x_data, y_data4_mean,'y',x_data, y_data4_mean,'*y');
% 
%     errorbar(x_data, y_data1_mean,y_data1_std);
%     plot(x_data, y_data1,'db');
%     
%     errorbar(x_data, y_data2_mean,y_data2_std);
%     plot(x_data, y_data2,'dr');
%     
%     errorbar(x_data, y_data3_mean,y_data3_std);
%     plot(x_data, y_data3,'dy');
%     
%     errorbar(x_data, y_data4_mean,y_data4_std);
%     plot(x_data, y_data4,'dm');
%     
%     grid minor;
%     legendsTemp = {'MOGA mean','MOGA max',...
%         'MOGA + NN mean','MOGA + NN mean max',...
%         'MOGA + re-scaling mean','MOGA + re-scaling max',...
%         'MOGA + NN + re-scaling mean','MOGA + NN + re-scaling max'};
%     legend(legendsTemp, 'Location', 'Southeast');
%     xlabel('Generation');   ylabel('fitness');
%     title(['Fitness "',fitnessOrder{1,i},'" over Generation']);
%     set(gca,'FontSize',12, 'FontWeight','bold');
%     hold off;
% end
end

%% Get parameters from the last generation and normalize them:
chosen_method = GA1;

paramRanges = MML.Gen.Range;
paramNames = MML.Gen.Keys;
allSeqs = chosen_method.GA.Seqs;
Param = zeros(size(allSeqs,1),size(allSeqs,2));
normParam = zeros(size(allSeqs,1),size(allSeqs,2));

for i=1:23
    Param(:,i) = allSeqs(:,i,last_gen);
    minParam = paramRanges(1,i);
    maxParam = paramRanges(2,i);
    normParam(:,i) = ( Param(:,i) - minParam ) ./ ( maxParam - minParam );
end

%% Choose two parameters to plot:
Param_names = seqOrder_extend;% fitnessOrder;
Param_name2Plot = {'tau','b'};

X1name = Param_name2Plot{1,1};
X2name = Param_name2Plot{1,2};
X1 = normParam(:,strcmp(X1name,Param_names));
X2 = normParam(:,strcmp(X2name,Param_names));

figure;
plot(X1,X2,'k*','MarkerSize',5);
title('param distribution');
xlabel(X1name);
ylabel(X2name);

%% MOOGA figures, show divercity in clustering:
% https://www.mathworks.com/help/stats/kmeans.html

num_of_clusters = 5;

opts = statset('Display','final');
[idx,C,sumd,D] = kmeans([X1,X2],num_of_clusters,'Distance','sqeuclidean',...
    'Replicates',5,'Options',opts);

%% plot clusters centers:
figure; hold on;
legends = cell(1,num_of_clusters*2);

for i=1:num_of_clusters
    plot(X1(idx==i,1),X2(idx==i,1),...
        colors{1,i},'MarkerSize',12)
    legends{1,2*i-1} = sprintf('Cluster %d',i);
    legends{1,2*i} = sprintf('Centroids %d',i);
    
    plot(C(i,1),C(i,2),'ko',...
     'MarkerSize',10,'LineWidth',4,...
     'MarkerFaceColor',colors1{1,i},'MarkerEdgeColor','k');
end
grid minor;
title(['clustering with ',num2str(num_of_clusters),' clusters']);
xlabel(X1name);
ylabel(X2name);
legend(legends{1,:});

%% Color the clusters areas:

%Use kmeans to compute the distance from each centroid to
% points on a grid. To do this, pass the centroids (C) and points on
% a grid to kmeans, and implement one iteration of the algorithm.

p1 = min(X1):0.001:max(X1);
p2 = min(X2):0.01:max(X2);
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
plot(X1,X2,'k*','MarkerSize',5);
xlabel(X1name);
ylabel(X2name);
legend(legends{1,1:2:end});
hold off;


%% display parameters divercity:
num_of_clusters = 5;

figure; hold on;
legends = cell(1,num_of_clusters*2);

chosen_method = {GA1,GA2,GA3,GA4};

Param_names = seqOrder_extend;
Param_name2Plot = {'tau','b'};
X1name = Param_name2Plot{1,1};
X2name = Param_name2Plot{1,2};

for j=1:4
    
    allSeqs = chosen_method{1,j}.GA.Seqs;
    Param = zeros(size(allSeqs,1),size(allSeqs,2));
    normParam = zeros(size(allSeqs,1),size(allSeqs,2));

    for i=1:23
        Param(:,i) = allSeqs(:,i,end);
        minParam = paramRanges(1,i);
        maxParam = paramRanges(2,i);
        normParam(:,i) = ( Param(:,i) - minParam ) ./ ( maxParam - minParam );
    end

    Xtemp = chosen_method{1,j}.GA.Seqs(:,:,end);

    X1 = normParam(:,strcmp(X1name,Param_names));
    X2 = normParam(:,strcmp(X2name,Param_names));

    opts = statset('Display','final');
    [idx,C,~,~] = kmeans([X1,X2],num_of_clusters,'Distance','sqeuclidean',...
        'Replicates',5,'Options',opts);
    
    Tend_ratio = chosen_method{1,j}.GA.Tend_ratio(:,1,last_gen);
    ids_good = (Tend_ratio >0.99);
    
    ax = subplot(2,2,j); hold on;
    for i=1:num_of_clusters
        plot(ax,X1(idx==i & ids_good,1),X2(idx==i & ids_good,1),...
            colors{1,i},'MarkerSize',18);
        plot(ax,X1(idx==i & ~ids_good,1),X2(idx==i & ~ids_good,1),...
            colors{1,i},'MarkerSize',5,'Marker','o','MarkerFaceColor','none');
        legends{1,2*i-1} = sprintf('Cluster %d',i);
        legends{1,2*i} = sprintf('Centroids %d',i);

        plot(ax,C(i,1),C(i,2),'ko',...
         'MarkerSize',10,'LineWidth',4,...
         'MarkerFaceColor',colors1{1,i},'MarkerEdgeColor','k');
    end
    
    grid minor;
    title(['clustering with ',num2str(num_of_clusters),...
        ' clusters: ',titleAdd{1,j}]);
    xlabel([X1name,' norm']);
    ylabel([X2name,' norm']);
    axis([0,1,0,1]);
    legend(legends{1,:},'Location','best');

end

%% display fitness divercity:
num_of_clusters = 5;

figure; hold on;
legends = cell(1,num_of_clusters*2);

chosen_method = {GA1,GA2,GA3,GA4};

method = 'Fit';
Param_names = fitnessOrder;
Param_name2Plot = {'VelFit','NrgEffFit'};

for j=1:4
    
    X = chosen_method{1,j}.GA.Fit(:,:,end);
    
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
        ' clusters: ',titleAdd{1,j}]);
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

inputs = ([Param1,Param2])';

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


%% display computation time:

chosen_method = {GA1,GA2,GA3,GA4};

totMLtime = zeros(4,10);
totSIMtime = zeros(4,10);

for j=1:4
    
    for k = 1:10
        allMLtimes = chosen_method{1,j}.GA.MLseqRunTime(:,1,k);
        allSIMtimes = chosen_method{1,j}.GA.simRunTime(:,1,k);

        totMLtime(j,k) = sum(allMLtimes,'omitnan');
        totSIMtime(j,k) = sum(allSIMtimes,'omitnan');
    end
    
end

MOOGAtime = totSIMtime(1,:);
MOOGA_NN_time = totSIMtime(2,:) + totMLtime(2,:);
MOOGA_rescale_time = totSIMtime(3,:) + 2 .* totMLtime(3,:); %because when rescale we run the sim twice
MOOGA_NN_rescale_time = totSIMtime(4,:) + 2 .* totMLtime(4,:);

figure; hold on;
plot(1:10,MOOGAtime);
plot(1:10,MOOGA_NN_time);
plot(1:10,MOOGA_rescale_time);
plot(1:10,MOOGA_NN_rescale_time);
grid minor;
title('total gen times over different methods');
xlabel('generation');
ylabel('time [sec]');
legend(legends,'Location','best');

% plot only the biPed sim time:
figure; hold on;
plot(1:10,totSIMtime(1,:));
plot(1:10,totSIMtime(2,:));
plot(1:10,totSIMtime(3,:));
plot(1:10,totSIMtime(4,:));
grid minor;
title('total biped sim time over generation num');
xlabel('generation');
ylabel('time [sec]');
legend(legends,'Location','best');

%% get sim endType:

fit_GA1 = GA1.GA.Fit(:,FitNum,:);

y_data1 = squeeze(max(fit_GA1,[],1));

%% try stuff:

X = GA1.GA.Tend_ratio(:,1,last_gen);
y_data1_max = squeeze(max(X,[],1));
y_data1_mean = squeeze(mean(X,1));

figure;
plot(1:25,y_data1_max); hold on;
plot(1:25,y_data1_mean);
errorbar(1:25,y_data1_mean,std(X,[],1));
legend('Max','mean');
title('MOGA only- no NN or rescale');
xlabel('generation num');
ylabel('T(end)/Sim.Tend');
grid minor


% X = GA1.GA.totGenTime;
% figure;
% plot(1:25,X);
% grid minor

