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
    % get y-axis data:
    y_data1 = squeeze(max(data{1,1}.GA.Fit(:,FitNum,1:gen_num),[],1));
    y_data2 = squeeze(max(data{1,2}.GA.Fit(:,FitNum,1:gen_num),[],1));
    y_data3 = squeeze(max(data{1,3}.GA.Fit(:,FitNum,1:gen_num),[],1));
    y_data4 = squeeze(max(data{1,4}.GA.Fit(:,FitNum,1:gen_num),[],1));

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

