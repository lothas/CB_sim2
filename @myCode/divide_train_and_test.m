function [obj] = divide_train_and_test(obj)
% this function divide the data to test group, validation group and train group

% ratio - the % of data in training group

switch size(obj.data_file_name,2)
    case 2 % divide the 1st file to train and valid and take the 2nd file as test
        samplesNum = obj.sampl_num_in_files(1,1);
        randIds = randsample(samplesNum,samplesNum);
        trainingSize = floor(obj.train_ratio*samplesNum);
        trainingIds = randIds(1:trainingSize);
        validationIds = randIds(trainingSize+1:end);
        
        obj.train_ind = trainingIds;
        obj.valid_ind = validationIds;
        obj.test_ind = 'the test group is the data from the 2nd file';
        
        [sampl_tr_val,targ_tr_val] = ...
            obj.prepareData_to_NN(obj.sim_results.results_train_and_valid,...
            obj.sim_periods.periods_train_and_valid,...
            obj.ids.ids_train_and_val);
        
        obj.sampl_train = sampl_tr_val(:,trainingIds);
        obj.targ_train = targ_tr_val(:,trainingIds);
        obj.sampl_valid = sampl_tr_val(:,validationIds);
        obj.targ_valid = targ_tr_val(:,validationIds);
        
        [sampl_test,targ_test] = ...
            obj.prepareData_to_NN(obj.sim_results.results_test,...
            obj.sim_periods.periods_test,...
            obj.ids.ids_test);
        obj.sampl_test = sampl_test;
        obj.targ_test = targ_test;
        
    case 3
        % then the 1st file will contain the train data and the 2nd one will
        % contain the validation data data and the 3rd the test data
        [sampl_tr,targ_tr] = ...
            obj.prepareData_to_NN(obj.sim_results.results_train,...
            obj.sim_periods.periods_train,...
            obj.ids.ids_train);
        obj.sampl_train = sampl_tr;
        obj.targ_train = targ_tr;
        
        [sampl_val,targ_val] = ...
            obj.prepareData_to_NN(obj.sim_results.results_valid,...
            obj.sim_periods.periods_valid,...
            obj.ids.ids_val);
        obj.sampl_valid = sampl_val;
        obj.targ_valid = targ_val;
        
        [sampl_test,targ_test] = ...
            obj.prepareData_to_NN(obj.sim_results.results_test,...
            obj.sim_periods.periods_test,...
            obj.ids.ids_test);
        obj.sampl_test = sampl_test;
        obj.targ_test = targ_test;
        
        obj.train_ind = 'the train group is the data from the 1st file';
        obj.valid_ind = 'the valid group is the data from the 2nd file';
        obj.test_ind = 'the test group is the data from the 3rd file';
        
    otherwise
        % is file name is a char, then size return the number...
        % of latter which is almost alway more than 3:)
        
        %divide the groups randomly
        % then divide the data randomly to train and test
        samplesNum = obj.sampl_num_in_files
        
        randIds = randsample(samplesNum,samplesNum);
        trainingSize = floor(obj.train_ratio*samplesNum);
        validSize = floor(obj.valid_ratio*samplesNum);

        trainingIds = randIds(1:trainingSize);
        validationIds = randIds((trainingSize+1):(validSize+trainingSize));
        testingIds = randIds(validSize+trainingSize:end);
        
        obj.train_ind = trainingIds;
        obj.valid_ind = validationIds;
        obj.test_ind = testingIds;
        
        [sampl,targ] = obj.prepareData_to_NN(obj.sim_results,obj.sim_periods,obj.ids);
        obj.sampl_train = sampl(:,trainingIds);
        obj.targ_train = targ(:,trainingIds);
        obj.sampl_valid = sampl(:,validationIds);
        obj.targ_valid = targ(:,validationIds);
        obj.sampl_test = sampl(:,testingIds);
        obj.targ_test = targ(:,testingIds);
end
  
% save the data in its undevided form for futuer bug checkig the the
% devision to train, validation and test.
switch size(obj.data_file_name,2)
    case 2
        obj.unGroupted_data.sampl_tr_val = sampl_tr_val;
        obj.unGroupted_data.targ_tr_val = targ_tr_val;
        obj.unGroupted_data.data_group_ind.train_ind = trainingIds;
        obj.unGroupted_data.data_group_ind.valid_ind = validationIds;
    case 3
        % there is no ungroup data (each group was divided manually).
    otherwise
        obj.unGroupted_data.sampl = sampl;
        obj.unGroupted_data.targ = targ;
        obj.unGroupted_data.data_group_ind.train_ind = trainingIds;
        obj.unGroupted_data.data_group_ind.valid_ind = validationIds;
        obj.unGroupted_data.data_group_ind.test_ind = testingIds;
end
        

end

