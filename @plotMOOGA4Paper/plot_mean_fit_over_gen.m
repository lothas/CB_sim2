function plot_mean_fit_over_gen(obj,whichFit,gen_num,whichGenes)
% this function plots the maximum fitness in every generation
%   and the Mean fitness (+error bars).

% 'whichFit' - row vector of the fitnesess to plot.
% 'gen_num' - up to which generation number to plot
% 'whichGenes' - which genes to take

x_data = 1:gen_num;
Markers = {'db','dr','dy','dm'};

for i = 1:length(whichFit)
    FitNum = whichFit(1,i);
    
    % open figure:
    figure; hold on;
    
    for j=1:4
        
        % extract the fitness:
        switch whichGenes
            case {'all','All'}
                y = obj.data{1,j}.GA.Fit(:,FitNum,x_data);
                
            case {'top_pop','topPop'}
                % take the top 15 percent
                popSize = obj.data{1,j}.GA.Population;
                topPopNum = floor(0.15 * popSize);
                [ID,~] = obj.data{1,j}.GA.GetTopPop(topPopNum);
                y = obj.data{1,j}.GA.Fit(ID,FitNum,x_data);
                
            otherwise
                error('invalid request');
        end
                
        % get y-axis data:
        switch obj.fitnessOrder{1,i}
            case {'VelFit','NrgEffFit','VelRangeFit #1','VelRangeFit #2',...
                    'VelRangeFit #4','VelRangeFit #6','EigenFit'}
                y_data = squeeze(max(y,[],1));
            case {'VelRangeFit #3','VelRangeFit #5',...
                    'VelRangeFit #7','VelRangeFit #8'} % because 's_slow" can be negative!
                y_data = squeeze(min(y,[],1));
            otherwise
                error('no such fitness');
        end
        
        % get Mean fit:
        y_data_mean = squeeze(mean(y,1));
        
        % get STdev fit:
        y_data_std = squeeze(std(y,[],1));
        
        % make plot:
        errorbar(x_data, y_data_mean,y_data_std);
        plot(x_data, y_data,Markers{1,j});
        
        
        
    end
    
    
    % make title::
        switch whichGenes
            case {'all','All'}
                Title = sprintf('%s over generations \n All samples',...
                    obj.fitnessOrder{1,FitNum});
            case {'top_pop','topPop'}
                Title = sprintf('%s over generations \n top 15%% samples',...
                    obj.fitnessOrder{1,FitNum});

        end
        
    grid minor;
    legendsTemp = {obj.legends{1,1},...
        [obj.legends{1,1},' best'],...
        obj.legends{1,2},...
        [obj.legends{1,2},' best'],...
        obj.legends{1,3},...
        [obj.legends{1,3},' best'],...
        obj.legends{1,4},...
        [obj.legends{1,4},' best']};
    
    legend(legendsTemp, 'Location', 'Southeast');
    xlabel('Generation');   ylabel('Fitness');
    title(Title);
    set(gca,'FontSize',12, 'FontWeight','bold');
    hold off;
end

end

