function plot_max_fit_over_gen(whichFit,fitnessOrder,...
    legends,data,gen_num)
% this function plots the maximum fitness in every generation

% 'whichFit' - row vector of the fitnesess to plot.
% 'fitnessOrder' - the order in which the fitnesses are orginaized
% 'legends' - the legend to plot from the graph
% 'data' - cell containing the GA data for the 4 cases
% 'gen_num' - up to which generation number to plot

x_data = 1:gen_num;


for i = 1:length(whichFit)
    FitNum = whichFit(1,i);
    
    y1 = data{1,1}.GA.Fit(:,FitNum,x_data);
    y2 = data{1,2}.GA.Fit(:,FitNum,x_data);
    y3 = data{1,3}.GA.Fit(:,FitNum,x_data);
    y4 = data{1,4}.GA.Fit(:,FitNum,x_data);
    
    % get y-axis data:
    switch fitnessOrder{1,i}
        case {'VelFit','NrgEffFit','VelRangeFit #1','VelRangeFit #2',...
                'VelRangeFit #4','VelRangeFit #5','VelRangeFit #6',...
                'VelRangeFit #7','VelRangeFit #8','EigenFit'}
            y_data1 = squeeze(max(y1,[],1));
            y_data2 = squeeze(max(y2,[],1));
            y_data3 = squeeze(max(y3,[],1));
            y_data4 = squeeze(max(y4,[],1));
        case {'VelRangeFit #3'} % because 's_slow" can be negative!
            y_data1 = squeeze(min(y1,[],1));
            y_data2 = squeeze(min(y2,[],1));
            y_data3 = squeeze(min(y3,[],1));
            y_data4 = squeeze(min(y4,[],1));
        otherwise
            error('no such fitness');
    end

    figure
    hold on
    h1=plot(x_data, y_data1);
    h2=plot(x_data, y_data2);
    h3=plot(x_data, y_data3);
    h4=plot(x_data, y_data4);
    grid minor;
    legend([h1,h2,h3,h4], legends, 'Location', 'Southeast');
    xlabel('Generation');   ylabel('Fitness');
    title(['Fitness "',fitnessOrder{1,i},'" over Generation']);
    set(gca,'FontSize',12, 'FontWeight','bold');
end

end

