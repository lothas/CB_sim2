function [obj] = trainNN(obj,showWindow)
% this function trains a neural network using NNtoolbox

% TODO: add more control on the training parameters

% 'showWindow' - wether to show the training GUI or not 

obj.NN.net.trainParam.showWindow = showWindow; % dont show training window

[obj.NN.net, obj.NN.net_perf] = train(obj.NN.net, obj.sampl, obj.targ);

end

