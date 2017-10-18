function [x_data,y_data_mean] = ...
    plot_mean_fit_over_gen(obj,whichFit,gen_num,whichGenes)
% this function plots the maximum fitness in every generation
%   and the Mean fitness (+error bars).

% 'whichFit' - row vector of the fitnesess to plot.
% 'gen_num' - up to which generation number to plot
% 'whichGenes' - which genes to take

x_data = 1:gen_num;
y_data_mean = zeros(numel(obj.data_names),gen_num);
y_data_std = zeros(numel(obj.data_names),gen_num);

for i = 1:length(whichFit)
    FitNum = whichFit(1,i);
    for j=1:numel(obj.data_names)
        
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
        
        % get Mean fit:
        y_data_mean(j,:) = squeeze(mean(y,1));
        
        % get STdev fit:
        y_data_std(j,:) = squeeze(std(y,[],1));
            
    end
    
    
    if nargout == 0   % make plots only the func is being called with no outputs
        % make title::
        switch whichGenes
            case {'all','All'}
                Title = sprintf('mean of %s over generations \n All samples',...
                    obj.fitnessOrder{1,FitNum});
            case {'top_pop','topPop'}
                Title = sprintf('mean of %s over generations \n top 15%% samples',...
                    obj.fitnessOrder{1,FitNum});

        end

        % % % open figure:
        figure; hold on;
        for j=1:numel(obj.data_names)
            plot(x_data, y_data_mean(j,:),'LineWidth',2);
            %errorbar(x_data, y_data_mean(j,:),y_data_std(j,:));
        end
        grid minor;
        legend(obj.Legends, 'Location', 'Best');
        xlabel('Generation');   ylabel('Fitness');
        title(Title);
        set(gca,'FontSize',12, 'FontWeight','bold');
        hold off;
    end
    
end

end

