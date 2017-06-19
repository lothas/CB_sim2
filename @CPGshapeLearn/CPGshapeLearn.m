classdef CPGshapeLearn
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        disp_information = true; %wether to display information or not
        
        % data from MatsuokaML sim
        data_file_name = [];
        sim_results = [];
        ids = []; % indecies of samples with not 'NaN' period and "good" data
        
        % names of inputs and outputs:
        inputsNames = [];
        outputNames = [];
        
        % NN groups devision
        sampl = [];
        targ = [];
        train_ind = []; % indecies for training group
        valid_ind = []; % indecies for validation group
        test_ind = []; % indecies for test group
        
        % NN data structure
        NN = []; % structure with fields:
                 %      net = 'feedfarward neural network;
                 %      net_perf = performance structure ('tr')
                 %      hiddenNeuronNum = hidden neurons vector
                 %      out_from_train = NN output from train group
                 %      out_from_test = NN output from test group
                 %      MSE_test_perf = NN perf (MSE) on test group
                 %      MSE_per_parameter = MSE performance for every
                 %          output element
    end
    
    methods
                 % %%% % Class constructor % %%% %
        function obj = CPGshapeLearn(varargin)
                    obj.data_file_name = varargin{1};
                    
                    fileName = varargin{1};
                    load(fileName,'results');
                    good_ids = obj.filter_bad_ids(results);
                    obj.sim_results = results;
                    obj.ids = good_ids;
                    
                    obj.inputsNames = varargin{2};
                    obj.outputNames = varargin{3};
                    
                   obj = obj.prepareData_NN();

                if nargin < 3
                    warning('invalid number of inputs');
                    error('input order should be: 1)name of data file 2)inputsNames 3)outputNames');
                end
        end
    
    end
    
end

