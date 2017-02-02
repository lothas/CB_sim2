function plotBars_with_errBars( means,stdevs,Names,label_Y,graph_title,graph_legend)
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
h=bar(means');
set(gca,'XTickLabel',Names);
set(h,'BarWidth',1);
for k=1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*k-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x, means(k,:), stdevs(k,:), 'k', 'linestyle', 'none');
end
legend(graph_legend);
ylabel(label_Y);
title(graph_title);
hold off

end

