function hist3D_CS(X1,X2,X3,Labels,Title,XYbinNum,zbinNum)
%HIST3D_CS 3D histogram in cross-sections
% 
% Inputs:
% *) 'X1' - vector containing the z-axes;
% *) 'X2','X3' - vectors containing the x and y axes
% *) 'Labels' - the names of the parameters
% *) 'XYbinNum','zbinNum' - number of bins in every direction
% 
% NOTE: the 'zbinNum' determine the number of cross-sections

% get the bin edges for the Z-axis:
[~,X1edges,X1bins] = histcounts(X1,zbinNum);

% define the bin edges for the X-Y axes:
X2edges = linspace(min(X2),max(X2),XYbinNum+1);
X3edges = linspace(min(X3),max(X3),XYbinNum+1);


Nums = zeros(XYbinNum,XYbinNum,zbinNum);
for i=1:zbinNum
    [Nums(:,:,i),~,~] =...
        histcounts2(X2(X1bins == i),X3(X1bins == i),X2edges,X3edges,...
        'Normalization','pdf');
end

% calculate bin centers:
X1center = X1edges(1:end-1) + diff(X1edges) / 2;
X2center = X2edges(1:end-1) + diff(X2edges) / 2;
X3center = X3edges(1:end-1) + diff(X3edges) / 2;

[X2g,X3g,X1g] = meshgrid(X2center,X3center,X1center);

slice(X2g,X3g,X1g,Nums,[],[],X1center);
ax= gca;

% find largest element in 'Nums':
Max = max(max(max(Nums)));

ax.CLim = [0,0.3*Max];
xlabel(Labels{1,2});
ylabel(Labels{1,3});
zlabel(Labels{1,1});
title(Title);

% binsNum = 5;
% [~,~,X1bins] = histcounts(X1,binsNum);
% for i=1:binsNum
%     figure;
%     histogram2(X2(X1bins == i),X3(X1bins == i),10,...
%         'DisplayStyle','tile','ShowEmptyBins','on',...
%         'Normalization','pdf');
%     title('2D distribution of "a" and "b"');
%     xlabel('a');
%     ylabel('b');
% end
end

