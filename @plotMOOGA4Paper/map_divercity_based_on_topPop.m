function map_divercity_based_on_topPop(obj,gen_num,num_of_clusters,...
    whichParam_from,type_from,whichParam_to,type_to)
% the the clusters Ids from Param space and map it to fitness space (and
% vice versa)
% 
% Inputs: 
% *) 'num_of_clusters' - the number of calusters for 'kmeans' function
% *) 'whichParam_from/to' - 1x2 cell array. can be anything from:
%      fitnessOrder = {'VelFit','NrgEffFit',...
%          'VelRangeFit #1','VelRangeFit #2','VelRangeFit #3',...
%          'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
%          'VelRangeFit #7','VelRangeFit #8','EigenFit'};
%      
%      seqOrder_extend = {'tau','b','c_1','c_2','c_3','c_4',...
%                       'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
%                       'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}',...
%                       'ks_\tau','ks_c1','ks_c2','ks_c3','ks_c4'};
% *) 'type_from/to' -either 'fit' for fitneses spcae,
%                     or 'param' for parametr space. 
% 
% 
% 


figure; hold on;

% Plot the original sapce:
for j=1:4
    X = [];  
    switch type_from
        case {'fit','Fit'}
            % than we are doing fitness divercity:
            % exstract the fitnesses:
            Y = obj.data{1,j}.GA.Fit(:,:,gen_num);
            for i=1:length(whichParam_from)
                X(:,i) = Y(:,strcmp(whichParam_from{1,i},obj.fitnessOrder));
                X_label{1,i} = whichParam_from{1,i};
            end
            X = X';
        case {'param','Param'}
            % than we are doint parameters divercity
            % normalize the seq parameters to be between '0' to '1':
            normParam = obj.param_norm_minMax(gen_num,j);

            % get the norm seq:
            for i=1:length(whichParam_from)
                X(i,:) =  normParam(:,strcmp(whichParam_from{1,i},...
                    obj.seqOrder_extend));
                X_label{1,i} = [whichParam_from{1,i},' norm'];
            end
            
        otherwise
            error('invalid type');
    end
   
    % % only fot the top 15% genes:
    howMany = floor(0.15 * size(X,2));
    [ID,~] = obj.data{1,j}.GA.GetTopPop(howMany);
    for i=1:length(whichParam_from)
        X_temp(i,:) = X(i,ID);
    end
    X = X_temp;

    % plot all that are left:
    ids_good = true(size(X,2),1);

    Title = sprintf('clustering with %d clusters:  %s \n for the Top 15%% of the population \n',...
                    num_of_clusters,obj.titleAdd{1,j});
    
    
    % clustering using 'kmeans':
    opts = statset('Display','final');
    [idx,C,~,~] = kmeans(X',num_of_clusters,'Distance','sqeuclidean',...
        'Replicates',5,'Options',opts);
    
    
    ax = subplot(2,2,j); hold on;
    
    idx_vec{1,j} = idx;
    
    Title = [Title,sprintf('num of points in each cluster = [')];
    for n=1:max(idx)
        Title = [Title,sprintf(' %d ',sum(idx==n))];
    end
    Title = [Title,sprintf(' ] \n Color order: [')];
    for n=1:max(idx)
        Title = [Title,sprintf(' %s , ',obj.colors{1,n})];
    end
    Title = [Title,sprintf(' ]')];
    
    obj.kMean_plot_clusters_and_centers(ax,X,...
        X_label,Title,C,idx,ids_good);
    
    axis([0,1,0,1]);
    

end

fileName = sprintf('the graph in which we did the clsutering');
savefig(fileName);

figure; hold on;


clear X_label
% map to the other space:
for j=1:4
    clear X Y ID X_temp
    switch type_to
        case {'fit','Fit'}
            % than we are doing fitness divercity:
            % exstract the fitnesses:
            Y = obj.data{1,j}.GA.Fit(:,:,gen_num);
            for i=1:length(whichParam_to)
                X(:,i) = Y(:,strcmp(whichParam_to{1,i},obj.fitnessOrder));
                X_label{1,i} = whichParam_to{1,i};
            end
            X = X';
        case {'param','Param'}
            % than we are doint parameters divercity
            % normalize the seq parameters to be between '0' to '1':
            normParam = obj.param_norm_minMax(gen_num,j);
            for i=1:length(whichParam_to)
                X(i,:) =  normParam(:,strcmp(whichParam_to{1,i},...
                    obj.seqOrder_extend));
                X_label{1,i} = [whichParam_to{1,i},' norm'];
            end
            
        otherwise
            error('invalid type');
    end
    
    % % only fot the top 15% genes:
    howMany = floor(0.15 * size(X,2));
    [ID,~] = obj.data{1,j}.GA.GetTopPop(howMany);
    for i=1:length(whichParam_to)
        X_temp(i,:) = X(i,ID);
    end
    X = X_temp;
    % plot all that are left:
    ids_good = true(size(X,2),1);

    Title = sprintf('map to fitness space with \n the same clusters IDs');
        
    
    ax = subplot(2,2,j); hold on;
    C = NaN(size(X));
    idx = idx_vec{1,j};
    obj.kMean_plot_clusters_and_centers(ax,X,X_label,Title,C,idx,ids_good);
    axis([0,1,0,1]);
    

end

fileName = sprintf('the graph in which we map the clusters ids');
savefig(fileName);
end

