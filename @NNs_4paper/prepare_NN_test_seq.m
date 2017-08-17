function seq_after_NN = prepare_NN_test_seq(obj,seq,inputsNames,targetsNames)
% this function prepare the matrices for the NN tesing in MOOGA

seq = seq';
seq_after_NN = seq;

sampl = zeros(length(inputsNames),size(seq,2));

% Prepare inputs:
for i = 1:length(inputsNames)
    p_name = inputsNames{1,i};
    switch p_name
        case 'period_desired'
            % randomly generate the desired periods
            per_des_Min = obj.MML.perLim(1,1);
            per_des_Max = obj.MML.perLim(1,2);
            sampl(i,:) = per_des_Min + ...
                ((per_des_Max-per_des_Min) * rand(1,length(seq)));
        case 'periods'
            error('cant be "periods", need to be "periods_desired"');
        otherwise
            sampl(i,:) = seq(strcmp(p_name,obj.seqOrder),:);
    end   
end

net = obj.NN.net;
theta_S1_new = net(sampl);

change_index = false(1,length(obj.seqOrder));
for i=1:length(targetsNames)
    change_index = change_index | strcmp(targetsNames{1,i},obj.seqOrder);
end

seq_after_NN(change_index,:) = theta_S1_new;

seq_after_NN = seq_after_NN';

end
