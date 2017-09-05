function [sampl,targ] = ...
    prepare_NN_data(inputsNames,outputsNames,...
    seqOrder,seq,periods)
%PREPARE_NN_DATA prepare NN inputs and outputs

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
            sampl(i,:) = seq(strcmp(p_name,seqOrder),:);
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
            targ(i,:) = seq(strcmp(p_name,seqOrder),:);
    end
end

end

