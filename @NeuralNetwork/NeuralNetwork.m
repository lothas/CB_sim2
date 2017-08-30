classdef NeuralNetwork < handle & matlab.mixin.Copyable
    % this class contains the Methods and properties of the NN tools to be
    % use with "NNs_4paper" class
    
    properties
        
        inputs = []; % NN inputs
        targets = []; % NN targets
        
        net = []; % the NN 
        net_train_function = 'trainlm'; % the NN training function
        architecture = []; % row vector with neural network architecture
        net_train_perf = []; % contain the structure of the NN training performancwe
        
        train_RMSE = []; % Root Mean Squar Error on train group
        valid_RMSE = []; % Root Mean Squar Error on valid group
        test_RMSE = []; % Root Mean Squar Error on test group
        
        train_R2 = []; % R^2 on train group
        valid_R2 = []; % R^2 on valid group
        test_R2 = []; % R^2 on on test group
        
        % slope of the regression graph on the different groups:
        train_slope = []; 
        valid_slope = []; 
        test_slope = []; 
        
        disp_information = true; % show garphs and stuff while training
    end
    
    methods
        function obj = NeuralNetwork(NNarchitecture,NNinputs,NNtargets)
            obj.architecture = NNarchitecture;
            obj.inputs = NNinputs;
            obj.targets = NNtargets;
        end
    end
    
end

