function [sampl,targ] = prepare_NN_train_data(obj,inputsNames,outputsNames)
% this function prepare the matrices for the NN training

seq = obj.seq(obj.osc_ids,:);
seq = seq';

periods = obj.periods(1,obj.osc_ids);

sampl = zeros(length(inputsNames),size(seq,2));
targ = zeros(length(outputsNames),size(seq,2));

% Prepare inputs:
for i = 1:length(inputsNames)
    p_name = inputsNames{1,i};
    switch p_name
        case {'periods','period'}
            sampl(i,:) = periods;
        case {'freq'}
            sampl(i,:) = 1./periods;
        otherwise
            sampl(i,:) = seq(strcmp(p_name,obj.seqOrder),:);
    end   
end
% Prepare inputs and outputs:
for i = 1:length(outputsNames)
    p_name = outputsNames{1,i};
    switch p_name
        case {'periods','period'}
            targ(i,:) = periods;
        case {'freq'}
            targ(i,:) = 1./periods;
        otherwise
            targ(i,:) = seq(strcmp(p_name,obj.seqOrder),:);
    end
end

end
