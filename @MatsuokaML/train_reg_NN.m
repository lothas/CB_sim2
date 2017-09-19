function [net, tr] = ...
    train_reg_NN(obj, samples, targets, architecture, NNSamples)
%TRAINNN Train a NN and calculate its performance
    netPerf = zeros(1, 4);    % Array to store NN performance
    sampPerf = zeros(1, NNSamples);
    sampPerfSc = zeros(1, NNSamples);

    % Create and train the NN
    net = feedforwardnet(architecture);
    [net, tr] = train(net, samples, targets);
    
end

