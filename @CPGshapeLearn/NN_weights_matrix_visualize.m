function NN_weights_matrix_visualize(obj,net,graphTitle)
% this function takes a NN and plot a color maps of its weights

% inputs:
% 1) 'net' - the NN structure
%
% outputs:
% none

parametersCells = obj.inputsNames;

weightsInput = net.IW;  
weightsInput = cell2mat(weightsInput);

weightsOutput = net.LW;
weightsOutput = cell2mat(weightsOutput);

bias = net.b;
bias = cell2mat(bias);

weights = horzcat(diag(weightsOutput)*weightsInput);  % multiply the output weights with the neurons outputs

[x,y] = meshgrid(1:size(weights,2),1:size(weights,1));   %# Create x and y coordinates for the strings
HiddenNumCells = num2cell((1:size(weights,1)));
textStrings = num2str(weights(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding

% figure;
imagesc(abs(weights));
colormap(gca,flipud(gray));  %# Change the colormap to gray (so higher values are
                         %#   black and lower values are white)
title(graphTitle);
% title(['Weights in a NN with ',num2str(length(HiddenNumCells)),' hidden neurons']);
hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(abs(weights(:)) > midValue,1,3);%# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
set(gca,'XTick',1:size(weights,2),...     %# Change the axes tick marks
        'XTickLabel',parametersCells,...  %#   and tick labels
        'YTick',1:size(weights,1),...
        'YTickLabel',HiddenNumCells,...
        'TickLength',[0 0]);

end

