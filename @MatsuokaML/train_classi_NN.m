function [net, tr] = ...
    train_classi_NN(obj, samples, targets, architecture)
%TRAIN_CLASSI_NN Train a classification NN 

    % NOTE: also optional to use autoEncoders toghether with "softmaxlayer"
    %   (it might work better)

    % Create and train the NN
    net = patternnet(architecture);
    [net, tr] = train(net, samples, targets);
             
    
end

