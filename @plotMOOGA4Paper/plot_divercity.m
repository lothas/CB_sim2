function clustersID = plot_divercity(obj,gen_num,num_of_clusters,...
    whichParam,normParam,plotType)
% this function plot the divercity on a 2D graph with 'kmeans'
% clusification
% 
% Inputs: 
% *) 'num_of_clusters' - the number of calusters for 'kmeans' function
% *) 'whichParam' - 1x2 cell array. can be anything from:
%      fitnessOrder = {'VelFit','NrgEffFit',...
%          'VelRangeFit #1','VelRangeFit #2','VelRangeFit #3',...
%          'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
%          'VelRangeFit #7','VelRangeFit #8','EigenFit'};
%      
%      seqOrder = {'tau','b',...};
% *) 'normParam' - 1Xlength(whichParam) boolean array. 
%                   '1'- norm;    '0'- don't norm
% *) 'plotType' - 
% 
% 
% 

% 1Xn cell array which stores the clusters indices for each case:
clustersID = cell(1,numel(obj.data_names));

% determine how many subplots:
    % if one GAfile     --> one plot
    % if between 2 to 4 --> 2X2
    % if between 5 to 9 --> 3X3         and ect...
numSP = ceil(sqrt(numel(obj.data_names)));

figure; hold on;

for j=1:numel(obj.data_names)
    X=[];
    X_label = whichParam;
    
    for p=1:length(whichParam)
        % check if which type of parameter is that:
        if any(strcmp(whichParam{1,p},obj.seqOrder))
            param_id = strcmp(whichParam{1,p},obj.seqOrder);
            P = obj.data{1,j}.GA.Seqs(:,param_id,gen_num);

            if normParam(1,p) % norm by genome min/Max
                minParam = obj.MML.Gen.Range(1,param_id);
                maxParam = obj.MML.Gen.Range(2,param_id);
                P = ( P - minParam ) ./ ( maxParam - minParam );
                X_label{1,p} = [whichParam{1,p},' norm'];
            end
            
        elseif any(strcmp(whichParam{1,p},obj.fitnessOrder))
            param_id = strcmp(whichParam{1,p},obj.fitnessOrder);
            P = obj.data{1,j}.GA.Fit(:,param_id,gen_num);
            
            if normParam(1,p) % norm by fitness min/Max
                P = ( P - min(P) ) ./ ( max(P) - min(P) );
                X_label{1,p} = [whichParam{1,p},' norm'];
            end
            
        end
        
        X(p,:) = P';
    end
    
    
    % waht type of plot we want:
    switch plotType
        case 'plot by Tend ratio'
            % get "good CPGs" (CBs that walking)by T(end)/Sim.Tend ratio
            Tend_ratio = obj.data{1,j}.GA.Tend_ratio(:,1,gen_num);
            ids_good = (Tend_ratio >0.95);
            
            Title = sprintf('clustering with %d clusters:  %s \n "o" for good CPG and "X" for a bad one \n based on Tend ratio',...
                        num_of_clusters,obj.Legends{1,j});
                    
        case 'plot by TopPop'
            % % only fot the top 15% genes:
            howMany = floor(0.15 * size(X,2));
            [ID,~] = obj.data{1,j}.GA.GetTopPop(howMany);
            for i=1:length(whichParam)
                X_temp(i,:) = X(i,ID);
            end
            X = X_temp;
            % plot all that are left:
            ids_good = true(size(X,2),1);
            
            Title = sprintf('clustering with %d clusters:  %s \n for the Top 15%% of the population',...
                            num_of_clusters,obj.Legends{1,j});
                        
        case 'plot by end condition'
            endCond = squeeze(obj.data{1,j}.GA.sim_endCond(:,1,gen_num));
            % '0' is the Sim.endCond that indicates that the sim reached to
            % t_span
            ids_good = (endCond == 0) | (endCond == 5) | (endCond == 6);
            
            Title = sprintf('clustering with %d clusters:  %s \n "o" Marker for good CPG and "X" for a bad one \n based on Sim.outType == 0,5,6',...
                        num_of_clusters,obj.Legends{1,j});
                    
        case 'plot by pareto fronts'
            % TODO: add this later
            
            
        otherwise
            error('invalid plotTypw');
    end
    
    % clustering using 'kmeans':
    opts = statset('Display','final');
    [idx,C,~,~] = kmeans(X',num_of_clusters,'Distance','sqeuclidean',...
        'Replicates',5,'Options',opts);
        
    ax = subplot(numSP,numSP,j); hold on;
    
    Title = [Title,sprintf('\n num of points in each cluster = [')];
    for n=1:max(idx)
        Title = [Title,sprintf(' %d ',sum(idx==n))];
    end
    Title = [Title,sprintf(' ]')];
    
    obj.kMean_plot_clusters_and_centers(ax,X,X_label,Title,C,idx,ids_good);
%     axis([0,1,0,1]);
    
    clustersID{1,j} = idx;
end

end

