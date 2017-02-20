function [obj] = shuffle_samples(obj)
% this function take the ungrouped data from "divide_train_and_test.m" and
% reshuffle the train/valid/test groups randomly (for cross validation
% purpocess)

% NOTE: this can only happen if there is a random devision to groups in the
% first place. meaning if we give only ONE data file to the class.

sampl = obj.unGroupted_data.sampl;
targ = obj.unGroupted_data.targ;

samplesNum = obj.sampl_num_in_files;

randIds = randsample(samplesNum,samplesNum);
trainingSize = floor(obj.train_ratio*samplesNum);
validSize = floor(obj.valid_ratio*samplesNum);

trainingIds = randIds(1:trainingSize);
validationIds = randIds((trainingSize+1):(validSize+trainingSize));
testingIds = randIds(validSize+trainingSize:end);

obj.sampl_train = sampl(:,trainingIds);
obj.targ_train = targ(:,trainingIds);
obj.sampl_valid = sampl(:,validationIds);
obj.targ_valid = targ(:,validationIds);
obj.sampl_test = sampl(:,testingIds);
obj.targ_test = targ(:,testingIds);

end

