function [ax_out] = heatmap_histogram(ax_in,vec,Edges,Title)
% this function get a vector (1 X N) and plot the histogram of that vector
% using a 1D heat map (a color bar)
% each color represent the height of the histogram bar

axes(ax_in);

[Prob,edges] = histcounts(vec,Edges, 'Normalization', 'probability');
bin_center = (edges(1:end-1)+ edges(2:end))/2;

h1 = imagesc(Prob);
caxis([min(Prob),max(Prob)]) % change color axis

% temp = h1.Parent.XTickLabel;
% h1.Parent.XTickLabel = {'0.7','0.8'};
set(h1.Parent,'XTick',bin_center,...
    'XTickLabel',num2str(bin_center));% something like this

h1.Parent.YTickLabel = '';

pbaspect([10 1 1]);
title(Title);

ax_out = gca;

end