function [sampl_train,targ_train,trainingIds,...
    sampl_test,targ_test,testingIds ] = divide_train_and_test(sampl,targ,ratio)
% this function divide the data to test group and train group

% ratio - the % of data in training group

randIds = randsample(size(sampl,2),size(sampl,2));
trainingSize = floor(ratio*size(sampl,2));
testSize = size(sampl,2) - trainingSize;
trainingIds = randIds(1:trainingSize);
testingIds = randIds(trainingSize:end);
sampl_train = sampl(:,trainingIds);
targ_train = targ(:,trainingIds);
sampl_test = sampl(:,testingIds);
targ_test = targ(:,testingIds);
end

