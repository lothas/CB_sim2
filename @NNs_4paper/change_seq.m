function [seq_after_NN,theta_S1_new] =...
    change_seq(obj,seq,method)
% this function use the NN/MoE to change the CPG sequece.

% create the inputs names for the test case
% (=using 'period_des' instead of 'period')
[inputsNames,targetsNames] = ...
    obj.check_NN_case(caseNum,'period_desired');

% save the Seq
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
        case 'period'
%             sampl(i,:) = periods;
              error('should not be periods');
        otherwise
            sampl(i,:) = seq(strcmp(p_name,obj.seqOrder),:);
    end   
end

theta_S1_new = obj.apply_net(sampl,method);

% check the the NN outputs are not crossing the allowed genome range:
if false % wether to limit the NN output
    for i=1:size(theta_S1_new,1)
        param_id = strcmp(targetsNames{1,i},obj.seqOrder);
        gen_min = obj.MML.Gen.Range(1,param_id);
        gen_max = obj.MML.Gen.Range(2,param_id);
        ids = (theta_S1_new(i,:) < gen_min) | (theta_S1_new(i,:) > gen_max);
        disp([num2str(sum(ids)),' are outside of the parameter range']);
        ind = find(ids);
        for j=1:sum(ids)
            theta_S1_new(i,ind(j)) =...
                min(max(theta_S1_new(i,ind(j)),gen_min),gen_max);
        end
    end
end

change_index = false(1,length(obj.seqOrder));
for i=1:length(targetsNames)
    change_index = change_index | strcmp(targetsNames{1,i},obj.seqOrder);
end

seq_after_NN(change_index,:) = theta_S1_new;

% Xedges = linspace(0,0.25,100);
% Yedges = linspace(0,10,100);
% figure;
% histogram2(theta_S1_new(1,:),theta_S1_new(2,:),Xedges,Yedges,...
%     'DisplayStyle','tile','ShowEmptyBins','on',...
%     'Normalization','pdf');
% title('2D distribution of tau and b after the NN');
% xlabel('tau');
% ylabel('b');
% axis([0,0.25,0,10]);

seq_after_NN = seq_after_NN';

end