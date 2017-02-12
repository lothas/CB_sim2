function [ExpertsWeights,gateWeights] = paper_MoE_init(obj,inputNum,outputCount,expertCount,dim)
% this function initialize new MoE network weights

% initialize parameters:
% initilize random weights by taking a random number from -1 to 1 devide by number of weights.
% gateNet weights-
gateWeights = (2*rand(expertCount, dim)-1)/inputNum; 
% experts weights-
ExpertsWeights = (2*rand(expertCount,dim, outputCount)-1)/inputNum;


end

