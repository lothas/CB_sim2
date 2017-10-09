function [ax] = plot_pointCloud_clusters(X,Y,Z,ids1,ids2,Names)
%PLOT_POINTCLOUD_CLUSTERS point each parametr on a 2D point cloud to try to
%identified the osc and n-osc clusters
% 
% Inputs:
% *) 'X','Y','Z' - three vectors
% *) 'ids1','ids2' - the ids of the clusters
% *) 'Names' - 1x3 cell array with the names of the three vectors


% make row vector with cluster ids:
clustes_ids = 1*ids1 + 2*ids2;

% % OPTION: don't show all samples
temp_id = randsample(length(X),50000);
X = X(1,temp_id);
Y = Y(1,temp_id);
Z = Z(1,temp_id);
clustes_ids = clustes_ids(1,temp_id);


% figure;
% gscatter(X,Y,{ids1,ids2});
% xlabel(Names{1,1});
% ylabel(Names{1,2});


% get the bin edges for the Z-axis:
zbinNum = 9;
[~,Zedges,Zbins] = histcounts(Z,zbinNum);

figure;
spNum = 1;  %subPlot number (starts with '1')
for i=1:max(Zbins)
     
    subplot(3,3,spNum);
    gscatter(X(1,Zbins==i),Y(1,Zbins==i),clustes_ids(1,Zbins==i));
    title({'cross-section of the parametrs point cloud at',...
        [num2str(Zedges(1,i)),' < ',Names{1,3},' < ',num2str(Zedges(1,i+1))]});
    xlabel(Names{1,1});
    ylabel(Names{1,2});
    legend('osc','n-osc');
    
    spNum = spNum + 1;
    
    % start a new subplot every 9 plots:
    if ~mod(i,9)
        spNum = 1;
        figure;
    end
end

ax = gca;

end

