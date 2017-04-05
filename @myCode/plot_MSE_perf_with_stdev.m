function plot_MSE_perf_with_stdev(obj,means,stdevs,Names,label_Y,graph_title,graph_legend)
% this function plots bars graphs with error bars centered on each bar.

% inputs:
% 1) means - a mean of a vector. the height of each bar
% 2) stdevs - standart deviation. size of error bar
% 3) Names - the 'x' axis labels
% 4) label_Y - 'y' axis label
% 5) graph_title - title of figure
% 6) graph_legend - legend

numgroups = size(means, 2); 
numbars = size(means, 1); 
groupwidth = min(0.8, numbars/(numbars+1.5));
figure; hold on

if numgroups == 1
    % if only one group of bars, make each bar different color
    h=bar(1, means(1),'b');
    hold on
    bar(2, means(2),'g');
else
    h=bar(means');
end

set(gca,'XTickLabel',Names);
set(h,'BarWidth',1);

for k=1:numbars
    if numgroups == 1
        errorbar(k, means(k,:), stdevs(k,:), 'k', 'linestyle', 'none');
    else
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*k-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x, means(k,:), stdevs(k,:), 'k', 'linestyle', 'none');
        % TODOL check if the errorbar is at the center of each bar.
    end
end

legend(graph_legend);
ylabel(label_Y);
title(graph_title);
hold off

end
