function plot_fit_over_gen(obj,whichFit,gen_num)
% this function plots the maximum fitness in every generation

% 'whichFit' - row vector of the fitnesess to plot.
% 'gen_num' - up to which generation number to plot

x_data = 1:gen_num;
y_data = zeros(4,gen_num);

for i = 1:length(whichFit)
    FitNum = whichFit(1,i);
    
    for j=1:numel(obj.data_names)
        % get the fit for all generations:
        y = obj.data{1,j}.GA.Fit(:,FitNum,x_data);
        
        % get wether this fit should be minimize/maximize:
            % '1'-fot max;      '-1'- for min
        minORmax = obj.data{1,j}.GA.FitMinMax;
        
        % get y-axis data:
        switch minORmax(1,j)
            case 1
                y_data(j,:) = squeeze(max(y,[],1));
            case -1
                y_data(j,:) = squeeze(min(y,[],1));
        end
        
    end
    
    figure; hold on;
    for j=1:numel(obj.data_names)
        plot(x_data, y_data(j,:));
    end
    grid minor;
    legend(obj.Legends, 'Location', 'Southeast');
    xlabel('Generation');   ylabel('Fitness');
    title(['Fitness "',obj.fitnessOrder{1,FitNum},'" over Generation']);
    set(gca,'FontSize',12, 'FontWeight','bold');

end
    
end

