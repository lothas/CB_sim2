function kMean_plot_clusters_and_centers(obj,ax,X1,X2,...
    X_label,Y_label,Title,centers,idx,goodSims_ids)
% takes results from 'kmeans' function and plot the clustering and 
%   the clusters' centes
% 
% Inputs:
% 'ax' - the axis to plot on.
% 'X1' 'X2' - the parameters to plot
% 'X_label' 'Y_label' 'Title' - axis labels and title
% 'centers' - the centers of the clusters.
% 'idx' - the classification indices for each data point to cluster
% 'goodSims_ids' - logic row vector containing the ids of the sims that
%                   "walked" all the way

num_of_clusters = size(centers,1);
ids = goodSims_ids;

% legends = cell(1,num_of_clusters*3);

for i=1:num_of_clusters
    % full Marker to "good" CPG:
    plot(ax,X1(idx==i & ids,1),X2(idx==i & ids,1),...
            obj.colors{1,i},'MarkerSize',18);
%     % empty Marker for "bad" CPG:
%     plot(ax,X1(idx==i & ~ids,1),X2(idx==i & ~ids,1),...
%         obj.colors{1,i},'MarkerSize',5,'Marker','o','MarkerFaceColor','none');
%     
    % 'X' Markers for "bad" CPG:
    plot(ax,X1(idx==i & ~ids,1),X2(idx==i & ~ids,1),...
        obj.colors{1,i},'MarkerSize',10,'Marker','X');
    
%     legends{1,2*i-1} = sprintf('Cluster %d',i);
%     legends{1,2*i} = sprintf('Centroids %d',i);
    
    plot(centers(i,1),centers(i,2),'ko',...
     'MarkerSize',10,'LineWidth',4,...
     'MarkerFaceColor',obj.colors1{1,i},'MarkerEdgeColor','k');
 
end

grid minor;
title(Title);
xlabel(X_label);
ylabel(Y_label);
% legend(legends{1,:},'Location','best');

% plot clusters areas:
% TODO: add it to independet function
if false
    %Use kmeans to compute the distance from each centroid to
    % points on a grid. To do this, pass the centroids (C) and points on
    % a grid to kmeans, and implement one iteration of the algorithm.

    p1 = min(X1):0.001:max(X1);
    p2 = min(X2):0.01:max(X2);
    [p1G,p2G] = meshgrid(p1,p2);
    P_Grid = [p1G(:),p2G(:)]; % Defines a fine grid on the plot

    idx2Region = kmeans(P_Grid,num_of_clusters,...
        'MaxIter',1,'Start',centers);
        % Assigns each node in the grid to the closest centroid

    % kmeans displays a warning stating that the algorithm did not
    % converge, which you should expect since the
    % software only implemented one iteration.    

    figure; hold on;
    for i=1:num_of_clusters
        plot(P_Grid(idx2Region==i,1),P_Grid(idx2Region==i,2),obj.colors{1,i});
    end
    plot(X1,X2,'k*','MarkerSize',5);
    xlabel(X1name);
    ylabel(X2name);
    legend(legends{1,1:2:end});
    hold off;
end

end

