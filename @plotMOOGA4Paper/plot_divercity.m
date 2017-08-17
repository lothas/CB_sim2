function plot_divercity(obj,gen_num,num_of_clusters,...
    whichParam,whichTypeOfParam,plotType)
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
%      seqOrder_extend = {'tau','b','c_1','c_2','c_3','c_4',...
%                       'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
%                       'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}',...
%                       'ks_\tau','ks_c1','ks_c2','ks_c3','ks_c4'};
% *) 'whichTypeOfParam' - can be either 'fitnessOrder' for fitness
%                           divercity. or 'seqOrder_extend' for parameters
%                           divercity.
% *) 'plotType' - 
% 
% 
% 

switch whichTypeOfParam % decide which sequqence to string compare later
    case {'fitness','fit','Fit'}
        order = obj.fitnessOrder;
    case {'param','Param'}
        order = obj.seqOrder_extend;
    otherwise
        error('invalid Type');
end

figure; hold on;

X1name = whichParam{1,1};
X2name = whichParam{1,2};

for j=1:4
    
%     allSeqs = data{1,j}.GA.Seqs;
    
    if any(strcmp('tau',order))
        % than we are doint parameters divercity
        % normalize the seq parameters to be between '0' to '1':
        normParam = obj.param_norm_minMax(gen_num,j);

        % get the norm seq:
        X1 = normParam(:,strcmp(X1name,order));
        X2 = normParam(:,strcmp(X2name,order));
        
        X_label = [X1name,' norm'];
        Y_label = [X2name,' norm'];
    else % than we are doing fitness divercity:
        % exstract the fitnesses:
        X = obj.data{1,j}.GA.Fit(:,:,gen_num);
        X1 = X(:,strcmp(X1name,order));
        X2 = X(:,strcmp(X2name,order));
        
        X_label = X1name;
        Y_label = X2name;
    end
    
    % waht type of plot we want:
    switch plotType
        case 'plot by Tend ratio'
            % get "good CPGs" (CBs that walking)by T(end)/Sim.Tend ratio
            Tend_ratio = obj.data{1,j}.GA.Tend_ratio(:,1,gen_num);
            ids_good = (Tend_ratio >0.95);
            
            Title = sprintf('clustering with %d clusters:  %s \n "o" for good CPG and "X" for a bad one \n based on Tend ratio',...
                        num_of_clusters,obj.titleAdd{1,j});
                    
        case 'plot by TopPop'
            % % only fot the top 15% genes:
            howMany = floor(0.15 * size(X1,1));
            [ID,~] = obj.data{1,j}.GA.GetTopPop(howMany);
            X1 = X1(ID,1);
            X2 = X2(ID,1);
            % plot all that are left:
            ids_good = true(length(X1),1);
            
            Title = sprintf('clustering with %d clusters:  %s \n for the Top 15%% of the population',...
                            num_of_clusters,obj.titleAdd{1,j});
                        
        case 'plot by end condition'
            endCond = squeeze(obj.data{1,j}.GA.sim_endCond(:,1,gen_num));
            % '0' is the Sim.endCond that indicates that the sim reached to
            % t_span
            ids_good = (endCond == 0) | (endCond == 5) | (endCond == 6);
            
            Title = sprintf('clustering with %d clusters:  %s \n "o" Marker for good CPG and "X" for a bad one \n based on Sim.outType == 0,5,6',...
                        num_of_clusters,obj.titleAdd{1,j});
                    
        case 'plot by pareto fronts'
            % TODO: add this later
            
            
        otherwise
            error('invalid plotTypw');
    end
    
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
end

