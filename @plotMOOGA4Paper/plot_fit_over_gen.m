function plot_fit_over_gen(obj,whichFit,gen_num)
% this function plots the maximum fitness in every generation

% 'whichFit' - row vector of the fitnesess to plot.
% 'gen_num' - up to which generation number to plot

x_data = 1:gen_num;
y_data = zeros(4,gen_num);

for i = 1:length(whichFit)
    FitNum = whichFit(1,i);
    
    for j=1:4
    y = obj.data{1,j}.GA.Fit(:,FitNum,x_data);
    
        % get y-axis data:
        switch obj.fitnessOrder{1,i}
            case {'VelFit','NrgEffFit','VelRangeFit #1','VelRangeFit #2',...
                    'VelRangeFit #4','VelRangeFit #6','EigenFit'}
                y_data(j,:) = squeeze(max(y,[],1));
            case {'VelRangeFit #3','VelRangeFit #5',...
                    'VelRangeFit #7','VelRangeFit #8'} % because 's_slow" can be negative!
                y_data(j,:) = squeeze(min(y,[],1));
            otherwise
                error('no such fitness');
        end
        
    end
    
    figure
    hold on
    h1=plot(x_data, y_data(1,:));
    h2=plot(x_data, y_data(2,:));
    h3=plot(x_data, y_data(3,:));
    h4=plot(x_data, y_data(4,:));
    grid minor;
    legend([h1,h2,h3,h4], obj.legends, 'Location', 'Southeast');
    xlabel('Generation');   ylabel('Fitness');
    title(['Fitness "',obj.fitnessOrder{1,i},'" over Generation']);
    set(gca,'FontSize',12, 'FontWeight','bold');

end
    
end

