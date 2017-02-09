function [obj] = divide_train_and_test(obj)
% this function divide the data to test group and train group

% ratio - the % of data in training group

samplesNum = size(obj.sampl,2);
randIds = randsample(samplesNum,samplesNum);
trainingSize = floor(obj.train2test_ratio*samplesNum);
% testSize = samplesNum - trainingSize;

trainingIds = randIds(1:trainingSize);
testingIds = randIds(trainingSize:end);

obj.train_ind = trainingIds;
obj.test_ind = testingIds;

obj.sampl_train = obj.sampl(:,trainingIds);
obj.targ_train = obj.targ(:,trainingIds);
obj.sampl_test = obj.sampl(:,testingIds);
obj.targ_test = obj.targ(:,testingIds);
end

