classdef myCode
    % this class contains all nof my functions for NN,MoE analisys
    
    properties
        
        % data from MatsuokaML sim
        data_file_name = [];
        Worig_or_What = false % if 'true' than training on the normalized weights
        results = []; 
        periods = [];
        ids = []; % indecies of samples with not 'NaN' period
        
        % names of inputs and outputs:
        seqOrder = []; % define the order in which the seq is coded
        inputsNames = [];
        outputNames = [];
        % 'b','tau' - beta and tau in the 4 neuron Matsuoka CPG (by defualt
        %               T=5*tau)
        % % only for 4 neuron CPG:
        %        'w_12',...,'w_43' - weights
        %       'W123',...,'W321' - combinations, for exampl 'W123' = 'w12'+'w23'+'w31' 
        %       'c_1',...'c_4' - tonic inputs
        %       and there are others like: 'sumW', 'prodW'
        % % only for 2 neuron symmetric CPG:
        % 'a','s' - are the weight and the tonic input
        % % possible outputs:
        % 'period','freq' (also every input can be defined as output).
        

        
        sampl = []; % the inputs to the NN
        targ = []; % the targets to the NN
        train2test_ratio = 0.8; % the ratio between training group and test group for the MoE training
        
        sampl_train = []; %training group
        targ_train = [];
        train_ind = []; % indecies for training group
        
        sampl_test = []; %test group
        targ_test = [];
        test_ind = []; % indecies for test group
        
        % parameters for use in "normal" nueral network:
        NN = []; % structure with fields:
                 %      net = 'feedfarward neural network;
                 %      net_perf = perfurmance structure ('tr')
                 %      hiddenNeuronNum = hidden neurons vector\
        
        % parameters for all MoE methods:
        expertCount = []; % the number of "experts" (each one is a NN)
        numOfIteretions = 100; % number of training cycles
        
        % parameters for use in our MoE method:
        my_MoE_out = []; % structure with fields:
                            %   'maxEphocs' - max number of ephocs for each NN training cycle
                            %   'ExpertHidLayer'- num of hidden layer in each expert
                            %   'ExpertHidNueron' - num of neurons in each hidden layer
                            %   'GateHidLayer' -num of hidden layer in gateNN
                            %   'GateHidNueron' - num of neurons in each hidden layer
                            %   'competetiveFlag' - our methods of training
                            %   'expertsNN' - contains all of the experts as NN structures
                            %   'gateNet' - contains the GateNN as a structure
                            %   'expertsTrainData' - structure contain the training data of the experts
                            %   'Moe_perf_over_iter' - overall MoE perf
                            %                           over interation num
                            %   'gateTraniData' - contain the training data of the gate
                            %   'out_from_train' - net outputs on train group
                            %   'out_from_test' - net outputs on test group
                            
                    % 'expertsTrainData' and 'gateTraniData' structure:
                            % 'expert_i_GroupSize' - the expert cluster size over interation num
                            % 'gateNN_perf_vec' - gateNet performance over iteration
                            % 'Experts_perf_mat' - experts best perf (MSE) over iteration num
                            % 'emptyGroupIndecator' - bolean matrix, '1'- every iteration that we have
                            %                           an empty group. for comparing the perf matrix
                            %                           sudden changes.

        % parameters for use in the paper MoE method:
        paper_MoE_out = []; % structure with fields:
                            %   'learningRate' - the learning rate 
                            %   'decay'- the decay of the learning rate 
        
    end
    
    methods
         % %%% % Class constructor % %%% %
        function obj = myCode(varargin)
            switch nargin
                case 4
                    obj.data_file_name = varargin{1};
                    
                    fileName = varargin{1};
                    obj = load_data(obj,fileName);
                    
                    obj.inputsNames = varargin{2};
                    obj.outputNames = varargin{3};
                    
                    obj.seqOrder =  varargin{4};
                    obj = prepareData_to_NN(obj);
                    obj = divide_train_and_test(obj);
                    
                otherwise
                    warning('invalid number of inputs');
                    error('input order should be: 1)name of data file 2)inputsNames 3)outputNames 4) order of seq');
            end
        end
        
        % %%% % load data to class % %%% %
        function [obj] = load_data(obj,fileName)
            load(fileName,'results');
            obj.results = results;
            periods = horzcat(results(:).periods);
            obj.periods = periods;
            % TODO: check if in the 4 neurons case the periods from results
            %        is matched.
            
            % get the Ids of only the ones with period (not 'NaN')
            ids_period = ~isnan(periods); 
            ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
            obj.ids = find(ids_period & ids_error);
        end
        
        % normalize the data by ( x_norm = (x-x_mean)/stdev(x) )  
        function [normalize_data,normParams] = normalizeData(data)
            normParams = zeros(size(data, 1), 2);
            normalize_data = zeros(size(data));
            for i = 1:size(data, 1)
                feat = data(i, :);
                normParams(i, :) = [mean(feat), std(feat)];
                normalize_data(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
            end
        end
    end
    
end

