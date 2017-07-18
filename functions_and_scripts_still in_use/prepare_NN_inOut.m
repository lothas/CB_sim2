function [sampl,targ] = prepare_NN_inOut(seq,periods,inputsNames,outputsNames,...
    seqOrder)
% this function prepare the matrices for the NN training

sampl = zeros(length(inputsNames),size(seq,2));
targ = zeros(length(outputsNames),size(seq,2));

% Prepare inputs:
for i = 1:length(inputsNames)
    p_name = inputsNames{1,i};
    switch p_name
        case 'period_desired'
            per_des_Min = 0.68;
            per_des_Max = 0.78;
            sampl(i,:) = per_des_Min + ...
                ((per_des_Max-per_des_Min) * rand(1,length(periods)));
            clear per_des_Min per_des_Max
        case 'periods'
            sampl(i,:) = periods;
        otherwise
            sampl(i,:) = seq(strcmp(p_name,seqOrder),:);
    end   
end
% Prepare inputs and outputs:
for i = 1:length(outputsNames)
    p_name = outputsNames{1,i};
    switch p_name
        case 'peridos'
            targ(i,:) = periods;
        otherwise
            targ(i,:) = seq(strcmp(p_name,seqOrder),:);
end

end

