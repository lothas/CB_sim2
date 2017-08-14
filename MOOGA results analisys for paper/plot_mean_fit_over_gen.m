function plot_mean_fit_over_gen(whichFit,fitnessOrder,...
    data,gen_num)
% this function plots the maximum fitness in every generation
%   and the Mean fitness (+error bars).

% 'whichFit' - row vector of the fitnesess to plot.
% 'fitnessOrder' - the order in which the fitnesses are orginaized
% 'data' - cell containing the GA data for the 4 cases
% 'gen_num' - up to which generation number to plot

x_data = 1:gen_num;


for i = 1:length(whichFit)
    FitNum = whichFit(1,i);
    
    fit_GA1 = data{1,1}.GA.Fit(:,FitNum,x_data);
    fit_GA2 = data{1,2}.GA.Fit(:,FitNum,x_data);
    fit_GA3 = data{1,3}.GA.Fit(:,FitNum,x_data);
    fit_GA4 = data{1,4}.GA.Fit(:,FitNum,x_data);
    
    % get y-axis data:
    y_data1 = squeeze(max(fit_GA1,[],1));
    y_data2 = squeeze(max(fit_GA2,[],1));
    y_data3 = squeeze(max(fit_GA3,[],1));
    y_data4 = squeeze(max(fit_GA4,[],1));

    % get Mean fit:
    y_data1_mean = squeeze(mean(fit_GA1,1));
    y_data2_mean = squeeze(mean(fit_GA2,1));
    y_data3_mean = squeeze(mean(fit_GA3,1));
    y_data4_mean = squeeze(mean(fit_GA4,1));

    % get STdev fit:
    y_data1_std = squeeze(std(fit_GA1,[],1));
    y_data2_std = squeeze(std(fit_GA2,[],1));
    y_data3_std = squeeze(std(fit_GA3,[],1));
    y_data4_std = squeeze(std(fit_GA4,[],1));

    figure; hold on;

    % plot GA only:
    errorbar(x_data, y_data1_mean,y_data1_std);
    plot(x_data, y_data1,'db');
    
    % plot GA + NN:
    errorbar(x_data, y_data2_mean,y_data2_std);
    plot(x_data, y_data2,'dr');
    
    % plot GA + rescale:
    errorbar(x_data, y_data3_mean,y_data3_std);
    plot(x_data, y_data3,'dy');
    
    % plot GA + NN + rescale:
    errorbar(x_data, y_data4_mean,y_data4_std);
    plot(x_data, y_data4,'dm');
    
    grid minor;
    legendsTemp = {'MOGA mean','MOGA max',...
        'MOGA + NN mean','MOGA + NN mean max',...
        'MOGA + re-scaling mean','MOGA + re-scaling max',...
        'MOGA + NN + re-scaling mean','MOGA + NN + re-scaling max'};
    legend(legendsTemp, 'Location', 'Southeast');
    xlabel('Generation');   ylabel('Fitness');
    title(['Fitness "',fitnessOrder{1,i},'" over Generation']);
    set(gca,'FontSize',12, 'FontWeight','bold');
    hold off;
end

end

