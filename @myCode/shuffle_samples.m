function [obj] = shuffle_samples(obj,mode)
% this function take the ungrouped data from "divide_train_and_test.m" and
% reshuffle the train/valid/test groups randomly (for cross validation
% purpocess)

% NOTE: this can only happen if there is a random devision to groups in the
% first place. meaning if we give only ONE data file to the class.

% sampl = obj.unGroupted_data.sampl;
% targ = obj.unGroupted_data.targ;
sampl = horzcat(obj.sampl_train,obj.sampl_valid,obj.sampl_test);
targ = horzcat(obj.targ_train,obj.targ_valid,obj.targ_test);

samplesNum = obj.sampl_num_in_files;

switch mode
    case 'completeShuffle'
        % shuffle all groups (traun,validation and test)
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
    case 'onlyTrainAndValid'
        % shuffle only the train and valiation groups
        trainSize = size(obj.sampl_train,2);
        validSize = size(obj.sampl_valid,2);
        trVal_group = trainSize+validSize;
        trainingIds = randsample(trVal_group,trainSize);
        validationIds = randsample(trVal_group,validSize);
        
        obj.sampl_train = sampl(:,trainingIds);
        obj.targ_train = targ(:,trainingIds);
        obj.sampl_valid = sampl(:,validationIds);
        obj.targ_valid = targ(:,validationIds);
        
    otherwise
        error('invalid mode');
end
        

end

