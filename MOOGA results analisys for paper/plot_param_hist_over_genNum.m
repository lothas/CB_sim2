function plot_param_hist_over_genNum(data,whichParam,gen_num,...
    fitnessOrder,seqOrder_extend,titleAdd)
% plot distribution of parameter over each generetion (to see hoe the GA is
% evolving)
% 
% Inputs: 
% *) 'data' - cell array contain the GA results
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
% *)'titleAdd' - for titles cases names 

X_data = 1:gen_num;
figure; 

for i=1:4
    switch whichParam
        case {'VelFit','NrgEffFit',...
             'VelRangeFit #1','VelRangeFit #2','VelRangeFit #3',...
             'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
             'VelRangeFit #7','VelRangeFit #8','EigenFit'}
            ID = strcmp(whichParam,fitnessOrder);
             Y = squeeze(data{1,i}.GA.Fit(:,ID,X_data));
            
            
        case {'tau','b','c_1','c_2','c_3','c_4',...
              'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
              'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}',...
              'ks_\tau','ks_c1','ks_c2','ks_c3','ks_c4'}
             ID = strcmp(whichParam,seqOrder_extend);
             Y = squeeze(data{1,i}.GA.Seqs(:,ID,X_data));
        otherwise
            error('invalid Param...');
    end
    
    absoluteMin = min(min(Y)); % to detemine the axis min
    absoluteMax = max(max(Y)); % to detemine the axis max
    
    lastBinEdge = absoluteMax*1.5;
    [Prob,edges] = distribution_over_genNum(Y,lastBinEdge);
    bin_center = (edges(1:end-1)+ edges(2:end))/2;
    
    subplot(2,2,i); hold on
    imagesc(X_data,bin_center,Prob');

    title({['the distribution of ',whichParam,' over genNum'],...
        ['for case: ',titleAdd{1,i}]});
    xlabel('generation num');
    ylabel(whichParam);
    grid minor
    axis([0,gen_num,absoluteMin*0.9,absoluteMax*1.1]);
    hold off

end

