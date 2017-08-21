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

X1name = whichParam_from{1,1};
X2name = whichParam_from{1,2};

% Plot the original sapce:
for j=1:4
      
    switch type_from
        case {'fit','Fit'}
            % than we are doing fitness divercity:
            % exstract the fitnesses:
            X = obj.data{1,j}.GA.Fit(:,:,gen_num);
            X1 = X(:,strcmp(X1name,obj.fitnessOrder));
            X2 = X(:,strcmp(X2name,obj.fitnessOrder));

            X_label = X1name;
            Y_label = X2name;
        case {'param','Param'}
            % than we are doint parameters divercity
            % normalize the seq parameters to be between '0' to '1':
            normParam = obj.param_norm_minMax(gen_num,j);

            % get the norm seq:
            X1 = normParam(:,strcmp(X1name,obj.seqOrder_extend));
            X2 = normParam(:,strcmp(X2name,obj.seqOrder_extend));

            X_label = [X1name,' norm'];
            Y_label = [X2name,' norm'];
        otherwise
            error('invalid type');
    end
   
    % % only fot the top 15% genes:
    howMany = floor(0.15 * size(X1,1));
    [ID,~] = obj.data{1,j}.GA.GetTopPop(howMany);
    X1 = X1(ID,1);
    X2 = X2(ID,1);
    % plot all that are left:
    ids_good = true(length(X1),1);

    Title = sprintf('clustering with %d clusters:  %s \n for the Top 15%% of the population',...
                    num_of_clusters,obj.titleAdd{1,j});
    
    
    % clustering using 'kmeans':
    opts = statset('Display','final');
    [idx,C,~,~] = kmeans([X1,X2],num_of_clusters,'Distance','sqeuclidean',...
        'Replicates',5,'Options',opts);
    
    
    ax = subplot(2,2,j); hold on;
    
    Title = [Title,sprintf('\n num of "good" points = %d',sum(ids_good))];
    obj.kMean_plot_clusters_and_centers(ax,X1,X2,...
        X_label,Y_label,Title,C,idx,ids_good);
    axis([0,1,0,1]);
    

end

figure; hold on;

X1name = whichParam_to{1,1};
X2name = whichParam_to{1,2};

% map to the other space:
for j=1:4
    
    switch type_to
        case {'fit','Fit'}
            % than we are doing fitness divercity:
            % exstract the fitnesses:
            X = obj.data{1,j}.GA.Fit(:,:,gen_num);
            X1 = X(:,strcmp(X1name,obj.fitnessOrder));
            X2 = X(:,strcmp(X2name,obj.fitnessOrder));

            X_label = X1name;
            Y_label = X2name;
        case {'param','Param'}
            % than we are doint parameters divercity
            % normalize the seq parameters to be between '0' to '1':
            normParam = obj.param_norm_minMax(gen_num,j);

            % get the norm seq:
            X1 = normParam(:,strcmp(X1name,obj.seqOrder_extend));
            X2 = normParam(:,strcmp(X2name,obj.seqOrder_extend));

            X_label = [X1name,' norm'];
            Y_label = [X2name,' norm'];
        otherwise
            error('invalid type');
    end
    
    % % only fot the top 15% genes:
    howMany = floor(0.15 * size(X1,1));
    [ID,~] = obj.data{1,j}.GA.GetTopPop(howMany);
    X1 = X1(ID,1);
    X2 = X2(ID,1);
    % plot all that are left:
    ids_good = true(length(X1),1);

    Title = sprintf('map to fitness space with \n the same clusters IDs');
        
    
    ax = subplot(2,2,j); hold on;
    C = NaN(num_of_clusters,2);
    obj.kMean_plot_clusters_and_centers(ax,X1,X2,...
        X_label,Y_label,Title,C,idx,ids_good);
    axis([0,1,0,1]);
    

end
end

