classdef MoE < handle & matlab.mixin.Copyable
    % "Mixture of Experts" regression model
    
    properties
        
        MoE_method = 'softCompetetetive';
        % MoE_method - *) 'hardCompetetetive' = "winner takes all"
%                      *) 'softCompetetetive' = "chance for everybody"
%                      *) 'collaboration' = (out = expertsOut*gateOut)
        
        expertCount = 3; % number of "experts" (NNs)
        architecture = 10; % row vector with neural network architecture
        inputs = []; % NN inputs
        targets = []; % NN targets
        
        train_Ind = []; % train group ids;
        valid_Ind = []; % validation group ids;
        test_Ind = []; % test group ids;
        
        numOfIteretions = 10; % number of training iterations
        MaxEpochs = 5; % maximum training epoch per itertation
        % each expert trains in total for "MaxEpochs" times the number of
        % iteration.
        
        expertsNN = []; % cell array cotains the experts as 'net' structure 
        gateNet = []; % gate NN (classifier)
        gateNN_perf_over_iter = [];
        
        train_MSE_vec = []; % MSE, train group, over iteration num
        valid_MSE_vec = []; % MSE, valid group, over iteration num
        test_MSE_vec = []; % MSE, test group, over iteration num
        
        train_RMSE = []; % Root Mean Squar Error on train group
        valid_RMSE = []; % Root Mean Squar Error on valid group
        test_RMSE = []; % Root Mean Squar Error on test group
        
        train_R2 = []; % Root Mean Squar Error on train group
        valid_R2 = []; % Root Mean Squar Error on valid group
        test_R2 = []; % Root Mean Squar Error on test group
        
        colors = []; % colors for the plotting functios
        legendNames = []; % for legend in the polotting functions
        
        % misc:
        gate_changes = []; %keep track on the gate chages per sample...
        
        disp_information = true; % show garphs and stuff while training
    end
    
    methods
        function obj = MoE(MoEinputs,MoEtargets,Method)
            obj.inputs = MoEinputs;
            obj.targets = MoEtargets;
            obj.MoE_method = Method;
            
            % TODO: need to devide manually for train/valid/test groups
            % otherwise, each iteration the groups whould be shuffled
        end
    end
    
end

