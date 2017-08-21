function kMean_plot_clusters_and_centers(obj,ax,X,...
    X_label,Title,centers,idx,goodSims_ids)
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



num_of_clusters = max(idx);
ids = goodSims_ids;

% legends = cell(1,num_of_clusters*3);
switch size(X,1)
    case 2 % 2D plot:
        for i=1:num_of_clusters
            % full Marker to "good" CPG:
            plot(ax,X(1,idx==i & ids),X(2,idx==i & ids),...
                    obj.colors{1,i},'MarkerSize',18);

            % 'X' Markers for "bad" CPG:
            plot(ax,X(1,idx==i & ~ids),X(2,idx==i & ~ids),...
                obj.colors{1,i},'MarkerSize',10,'Marker','X');
            
            if any(~isnan(centers))
                plot(centers(i,1),centers(i,2),'ko',...
                 'MarkerSize',10,'LineWidth',4,...
                 'MarkerFaceColor',obj.colors1{1,i},'MarkerEdgeColor','k');
            end
        end
        grid minor;
        title(Title);
        xlabel(X_label{1,1});
        ylabel(X_label{1,2});
        
    case 3 % 3D scatter
        for i=1:num_of_clusters
            X1_good = X(1,idx==i & ids);
            X2_good = X(2,idx==i & ids);
            X3_good = X(3,idx==i & ids);
            
            X1_bad = X(1,idx==i & ~ids);
            X2_bad = X(2,idx==i & ~ids);
            X3_bad = X(3,idx==i & ~ids);
        % full Marker to "good" CPG:
            scatter3(ax,X1_good,X2_good,X3_good,40,obj.colors{1,i},...
                'Marker','o');

            % 'X' Markers for "bad" CPG:
            scatter3(ax,X1_bad,X2_bad,X3_bad,40,obj.colors{1,i},...
                'Marker','x');
            
%             scatter3(X1_good, X2_good, 0*X3_good);    %projection from Z+
%             scatter3(X1_good, 0*X2_good, X3_good);    %projection from Y-
%             scatter3(0*X1_good, X2_good, X3_good);    %projection from X+

            if any(~isnan(centers))
                scatter3(centers(i,1),centers(i,2),centers(i,3),100,'ko',...
                 'LineWidth',4,...
                 'MarkerFaceColor',obj.colors1{1,i},'MarkerEdgeColor','k');
            end
        end
        grid minor;
        
        title(Title);
        xlabel(X_label{1,1});
        ylabel(X_label{1,2});
        zlabel(X_label{1,3});
        
    otherwise
        error('invalid number of parameters')
end          

end

