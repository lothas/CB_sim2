function randOut = rand_from_hist(N,histEdges,binHeight)
%RAND_FROM_HIST ruffle a random number given a specific histogram

% Inputs:
% *) 'N' - the number of points to get
% *) 'histEdges' - the edges given from the histogram 
% *) 'binHeight' - the height of the bin (units of probability)
% NOTE: you must use: "histcounts(X,'Normalization','probability')"

% comulative distribution function:
CDF = cumsum(binHeight);

% get N random numbers in a uniform range between '0' to '1':
R = rand(N,1);

% define a function for finding the probability from the CDF:
p = @(R) find(R < CDF, 1, 'first');

% get the binNum from which each 'R' belongs to:
%   the higher the probability (binHeight) the higher the change for the
%   bin to get selected.
rR = arrayfun(p,R);

% define a function to ruffle a random number given the bin num
%   the function finds the bin edges and ruffle the randNum between the
%   edges
get_from_bin = @(ind) (histEdges(ind-1) + rand()*...
    (histEdges(ind) - histEdges(ind-1)));

randOut= arrayfun(get_from_bin,rR);

end

