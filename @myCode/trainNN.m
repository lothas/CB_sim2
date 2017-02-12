function [obj] = trainNN(obj,showWindow)
% this function trains a neural network using NNtoolbox

% TODO: add more control on the training parameters

% 'showWindow' - wether to show the training GUI or not 

obj.NN.net.trainParam.showWindow = showWindow; % dont show training window

[obj.NN.net, obj.NN.net_perf] = train(obj.NN.net, obj.sampl_train, obj.targ_train);

obj.NN.out_from_train = obj.NN.net(obj.sampl_train);
obj.NN.out_from_test = obj.NN.net(obj.sampl_test);
end

