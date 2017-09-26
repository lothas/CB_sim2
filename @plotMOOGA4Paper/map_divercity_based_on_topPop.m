function map_divercity_based_on_topPop(obj,gen_num,num_of_clusters,...
    whichParam_from,normParam_from,whichParam_to,normParam_to)
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


% determine how many subplots:
    % if one GAfile     --> one plot
    % if between 2 to 4 --> 2X2
    % if between 5 to 9 --> 3X3         and ect...
numSP = ceil(sqrt(numel(obj.data_names)));

% Plot the original sapce:
clustersID = obj.plot_divercity(gen_num,num_of_clusters,...
    whichParam_from,normParam_from,'plot by TopPop');

% fileName = sprintf('the graph in which we did the clsutering');
% savefig(fileName);

figure; hold on;


clear X_label
% map to the other space:
for j=1:numel(obj.data_names)
    clear X ID X_temp
    
    X=[];
    X_label = whichParam_to;
    
    for p=1:length(whichParam_to)
        % check if which type of parameter is that:
        if any(strcmp(whichParam_to{1,p},obj.seqOrder))
            param_id = strcmp(whichParam_to{1,p},obj.seqOrder);
            P = obj.data{1,j}.GA.Seqs(:,param_id,gen_num);

            if normParam_to(1,p) % norm by genome min/Max
                minParam = obj.MML.Gen.Range(1,param_id);
                maxParam = obj.MML.Gen.Range(2,param_id);
                P = ( P - minParam ) ./ ( maxParam - minParam );
                X_label{1,p} = [whichParam_to{1,p},' norm'];
            end
            
        elseif any(strcmp(whichParam_to{1,p},obj.fitnessOrder))
            param_id = strcmp(whichParam_to{1,p},obj.fitnessOrder);
            P = obj.data{1,j}.GA.Fit(:,param_id,gen_num);
            
            if normParam_to(1,p) % norm by fitness min/Max
                P = ( P - min(P) ) ./ ( max(P) - min(P) );
                X_label{1,p} = [whichParam_to{1,p},' norm'];
            end
            
        end
        
        X(p,:) = P';
    end
    
    % % only for the top 15% genes:
    howMany = floor(0.15 * size(X,2));
    [ID,~] = obj.data{1,j}.GA.GetTopPop(howMany);
    for i=1:length(whichParam_to)
        X_temp(i,:) = X(i,ID);
    end
    X = X_temp;
    % plot all that are left:
    ids_good = true(size(X,2),1);

    Title = sprintf('map to fitness space with \n the same clusters IDs');
        
    
    ax = subplot(numSP,numSP,j); hold on;
    C = NaN(size(X));
    idx = clustersID{1,j};
    obj.kMean_plot_clusters_and_centers(ax,X,X_label,Title,C,idx,ids_good);
%     axis([0,1,0,1]);
    

end

% fileName = sprintf('the graph in which we map the clusters ids');
% savefig(fileName);
end

