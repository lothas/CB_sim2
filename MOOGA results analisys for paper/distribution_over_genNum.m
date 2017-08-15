function [Prob,edges] = distribution_over_genNum(Y,lastBinEdge)
% this function prepare a matrix with 'm' number of histigrams of a given
% parameter (one for each GA generation).
% 
% Inputs:
% *) 'Y' - (n x m) matrix
%           'n' - number of samples in each generation
%           'm' - number of generation
% 
% Outputs:
% *) 'Prob' - (n x m) matrix with the probability historgram
% 

N = 100; % the number of bins
Edges = linspace(0,lastBinEdge,N); % the bin edges

genNum = size(Y,2); % generations number

Prob = zeros(genNum,N-1);
for j=1:genNum
    Y_curr = Y(:,j);
    [Prob(j,:),edges] = histcounts(Y_curr,Edges,...
        'Normalization', 'probability');
    % bin_center = (edges(1:end-1)+ edges(2:end))/2;
end 


end

