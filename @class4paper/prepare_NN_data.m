function [sampl,targ] = prepare_NN_data(obj,ids,...
    inputsNames,outputsNames)
% this function prepare the matrices for the NN training

seq = obj.results(ids);
periods = mean(horzcat(results(ids).periods),1);

sampl = zeros(length(inputsNames),size(seq,2));
targ = zeros(length(outputsNames),size(seq,2));

% Prepare inputs:
for i = 1:length(inputsNames)
    p_name = inputsNames{1,i};
    switch p_name
        case 'period_desired'
            % randomly generate the desired periods
            per_des_Min = obj.MML.perLim(1,1);
            per_des_Max = obj.MML.perLim(1,2);
            sampl(i,:) = per_des_Min + ...
                ((per_des_Max-per_des_Min) * rand(1,length(periods)));
            clear per_des_Min per_des_Max
        case {'periods','period'}
            sampl(i,:) = periods;
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
        otherwise
            targ(i,:) = seq(strcmp(p_name,obj.seqOrder),:);
end

end
