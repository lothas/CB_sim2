function plot_param_hist_over_genNum(obj,whichParam,gen_num)
% plot distribution of parameter over each generetion (to see hoe the GA is
% evolving)
% 
% Inputs: 
% *) 'whichParam' - can be anything from:
%      fitnessOrder = {'VelFit','NrgEffFit',...
%          'VelRangeFit #1','VelRangeFit #2','VelRangeFit #3',...
%          'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
%          'VelRangeFit #7','VelRangeFit #8','EigenFit'};
%      
%      seqOrder_extend = {'tau','b','c_1','c_2','c_3','c_4',...
%                       'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
%                       'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}',...
%                       'ks_\tau','ks_c1','ks_c2','ks_c3','ks_c4'};
% *) 'gen_num' - the final generation that we want to plot     

% determine how many subplots:
    % if one GAfile     --> one plot
    % if between 2 to 4 --> 2X2
    % if between 5 to 9 --> 3X3         and ect...
numSP = ceil(sqrt(numel(obj.data_names)));

X_data = 1:gen_num;
figure; 

for i=1:numel(obj.data_names)
    % check if which type of parameter is that:
    if any(strcmp(whichParam,obj.seqOrder))
        param_id = strcmp(whichParam,obj.seqOrder);
        Y = squeeze(obj.data{1,i}.GA.Seqs(:,param_id,X_data));

    elseif any(strcmp(whichParam,obj.fitnessOrder))
        param_id = strcmp(whichParam,obj.fitnessOrder);
        Y = squeeze(obj.data{1,i}.GA.Fit(:,param_id,X_data));
    end
    
    absoluteMin = min(min(Y)); % to detemine the axis min
    absoluteMax = max(max(Y)); % to detemine the axis max
    
    lastBinEdge = absoluteMax*1.5;
    [Prob,edges] = obj.distribution_over_genNum(Y,lastBinEdge);
    bin_center = (edges(1:end-1)+ edges(2:end))/2;
    
    subplot(numSP,numSP,i); hold on
    
    imagesc(X_data,bin_center,Prob');

    title({['the distribution of "',whichParam,'" over genNum'],...
        ['for case: ',obj.Legends{1,i}]});
    xlabel('generation num');
    ylabel(whichParam);
    grid minor
    axis([0,gen_num,absoluteMin*0.9,absoluteMax*1.1]);
    hold off

end

