function [obj] = divide_train_and_test(obj)
% this function divide the data to test group and train group

% ratio - the % of data in training group

samplesNum = size(obj.sampl,2);

if size(obj.data_file_name,2) == 1
    % then divide the data randomly to train and test
    randIds = randsample(samplesNum,samplesNum);
    trainingSize = floor(obj.train2test_ratio*samplesNum);
    % testSize = samplesNum - trainingSize;

    trainingIds = randIds(1:trainingSize);
    testingIds = randIds(trainingSize:end);

else
    % then the 1st file will contain the train data and the 2nd one will
    % contain the test data
    trainingIds = 1:obj.sampl_num_in_files(1,1);
    testingIds = obj.sampl_num_in_files(1,1)+2:...
        (obj.sampl_num_in_files(1,1)+obj.sampl_num_in_files(1,2));
end
    
obj.train_ind = trainingIds;
obj.test_ind = testingIds;
    
obj.sampl_train = obj.sampl(:,trainingIds);
obj.targ_train = obj.targ(:,trainingIds);
obj.sampl_test = obj.sampl(:,testingIds);
obj.targ_test = obj.targ(:,testingIds);
end

