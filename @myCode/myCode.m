classdef myCode < handle & matlab.mixin.Copyable
    % this class contains all nof my functions for NN,MoE analisys
    
    properties
        
        sizeOfCPG=[]; % the number of neurons in the Matsuoka CPG
        rearrenge_in_a_uniq_way = [] % decide wether is needed to rearenge the genes in a uniq way for training
        
        disp_information = true; %wether to display information or not
        
        % data from MatsuokaML sim
        data_file_name = [];
        sampl_num_in_files = []; % if more than one file. than contain the number of samples in each file
        Worig_or_What = false % if 'true' than training on the normalized weights
        sim_results = [];   
        sim_periods = [];
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
        

        unGroupted_data = [] % save the data in its undevided form for futuer bug checkig the the
                             % devision to train, validation and test.
        train_ratio = 0.75; % the % of samples in the train group
        valid_ratio = 0.15; % the % of samples in the validation group
        %note: the rest will go to the test group
        
        sampl_train = []; %training group
        targ_train = [];
        train_ind = []; % indecies for training group
        
        sampl_valid = []; %validation group
        targ_valid = [];
        valid_ind = []; % indecies for validation group
        
        sampl_test = []; %test group
        targ_test = [];
        test_ind = []; % indecies for test group
        
        % parameters for use in "normal" nueral network:
        NN = []; % structure with fields:
                 %      net = 'feedfarward neural network;
                 %      net_perf = perfurmance structure ('tr')
                 %      hiddenNeuronNum = hidden neurons vector
                 %      out_from_train = NN output from train group
                 %      out_from_test = NN output from test group
                 %      MSE_test_perf = NN perf (MSE) on test group
        
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
        
        colors = []; % the colors of each expert on the graphs
        legendNames = []; % legend string array
    end
    
    methods
         % %%% % Class constructor % %%% %
        function obj = myCode(varargin)
                    obj.data_file_name = varargin{1};
                    
                    fileName = varargin{1};
                    obj = load_data(obj,fileName);
                    
                    obj.inputsNames = varargin{2};
                    obj.outputNames = varargin{3};
                    
                    obj.seqOrder =  varargin{4};
                    
                    obj.sizeOfCPG = varargin{5};
            switch nargin
                case 5
                    obj.rearrenge_in_a_uniq_way = false; 
                case 6
                    
                    obj.rearrenge_in_a_uniq_way = varargin{6};
                otherwise
                    warning('invalid number of inputs');
                    error('input order should be: 1)name of data file 2)inputsNames 3)outputNames 4) order of seq');
            end
            obj = divide_train_and_test(obj);
        end
        
        % %%% % load data to class % %%% %
        function [obj] = load_data(obj,fileName) 
            switch size(fileName,2)
                case 2
                    % if we give 2 file names then the 1st one will
                    %devide randomly to train and validation
                    % and the 2nd will be the test.
                    
                    % prepare the vector with the the amount of "good ids"
                    % in each data file:
                    sampl_num_in_files_temp = zeros(1,2);
                    
                    load(fileName{1,1},'results')
                    obj.sim_results.results_train_and_valid = results;
                    [periods_tr_val,ids_tr_val] = obj.filter_Nan_periods_and_get_ids(results);
                    obj.sim_periods.periods_train_and_valid = periods_tr_val;
                    obj.ids.ids_train_and_val = ids_tr_val;
                    sampl_num_in_files_temp(1,1) = length(ids_tr_val);
                    clear results
                    
                    load(fileName{1,2},'results')
                    obj.sim_results.results_test = results;
                    [periods_test,ids_test] = obj.filter_Nan_periods_and_get_ids(results);
                    obj.sim_periods.periods_test = periods_test;
                    obj.ids.ids_test = ids_test;
                    sampl_num_in_files_temp(1,2) = length(ids_test);
                    clear results
                                   
                case 3
                    % if we give 3 file names then the 1st one will be the
                    % train, the 2nd one will be validation and the 3rd
                    % will be the test
                    sampl_num_in_files_temp = zeros(1,3);
                    
                    load(fileName{1,1},'results')
                    obj.sim_results.results_train = results;
                    [periods_tr,ids_tr] = obj.filter_Nan_periods_and_get_ids(results);
                    obj.sim_periods.periods_train = periods_tr;
                    obj.ids.ids_train = ids_tr;
                    sampl_num_in_files_temp(1,1) = length(ids_tr);
                    clear results
                    
                    load(fileName{1,2},'results')
                    obj.sim_results.results__valid = results;
                    [periods_val,ids_val] = obj.filter_Nan_periods_and_get_ids(results);
                    obj.sim_periods.periods_valid = periods_val;
                    obj.ids.ids_val = ids_val;
                    sampl_num_in_files_temp(1,2) = length(ids_val);
                    clear results
                    
                    load(fileName{1,3},'results')
                    obj.sim_results.results_test = results;
                    [periods_test,ids_test] = obj.filter_Nan_periods_and_get_ids(results);
                    obj.sim_periods.periods_test = periods_test;
                    obj.ids.ids_test = ids_test;
                    sampl_num_in_files_temp(1,3) = length(ids_test);
                    clear results
                otherwise % is file name is a char, then size return the number of latter which is almost alway more than 3:)
                    load(fileName,'results');
                    obj.sampl_num_in_files = length(results);
                    [periods,good_ids] = obj.filter_Nan_periods_and_get_ids(results);
                    obj.sim_periods = periods;
                    obj.sim_results = results;
                    obj.ids = good_ids;
                    
                    sampl_num_in_files_temp = length(good_ids);
                    
            end
            obj.sampl_num_in_files = sampl_num_in_files_temp;
        end

        
        function [periods,ids] = filter_Nan_periods_and_get_ids(obj,results)
            % extract the periods:
            periods = horzcat(results(:).periods);
            
            % TODO: check if in the 4 neurons case the periods from results
            %        is matched.
            
            % get the Ids of only the ones with period (not 'NaN')
            ids_period = ~isnan(periods);
            % only ones with low enought error
            ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)';
            
            % for the 4Nuerons case the period contain two rows,
            %    one for the hip signal and one for the ankle signal:
            if (size(periods,1)>1) 
                % take only sample with the same periods in ankle and hip
                ids_similiar_periods = (periods(1,:)-periods(2,:)) < 0.1;
                % take only samples that don't have NaN in any joint
                ids_period = ids_period(1,:) & ids_period(2,:);
                % combain the two:
                ids_period = ids_similiar_periods & ids_period;
                
            end
            ids = find(ids_period & ids_error);         
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

