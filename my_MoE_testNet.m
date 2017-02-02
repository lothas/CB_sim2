function [netOut,gateOut,targ,belongToExpert,cluster_i_train_ind] = my_MoE_testNet(NNinputs,NNtargets,expertsNN,...
    gateNet,competetiveFlag)
%this function train a "Mixture of Experts" newural networks

%%%%%%%%% TODO: make this code to work with multiple NN outputs! %%%%%%%

% inputs:
% 1) NNinputs - the inputs for the neural networks,
%               matrix of (#inputs x #samples) dimention
% 2) NNtargets - the tarhets for the neural networks,
%                matrix of (#targets x #samples) dimention
% 3) expertsNN - cell array containing the experts and their performance.
% 4) gateNet - the Gate NN.
% 5) competetiveFlag - if '1'- "winner takes all"
%                         '2'- "chance for everybody"
%                         '3'- out = expertsOut * gateOut

% outputs:
% 1) netOut - output from MoE
% 2) gateOut - the gate output, represents the probability of a sample to
%               belong to each expert
% 3) targ - targets, rearenge the vector just in case of data orginizing
% 4) belongToExpert - which expert's cluster the sample belongs.
% 5) cluster_i_train_ind - in case of clustering. this is the indecies
%   specify which samples went to which expert.

expertCount = size(expertsNN,2);

belongToExpert = zeros(1,size(NNinputs,2));
expertsOut = zeros(expertCount,size(NNinputs,2)); % Experts output matrix 
netOut = [];
targ = [];
cluster_i_train_ind = cell(1,expertCount);

gateOut = gateNet(NNinputs);

switch competetiveFlag
    case 1
       % each data sample is classified based on the highest probability. 
       % for expample, given 2 experts and a datasample with
       % p1=0.6 chance to belong to the 1st expert and p2=0.4 chance to
       % belong to the 2nd expert. so this sample belongs to 1st expert. 
       [~,belongToExpert] = max(gateOut,[],1); % indecies of points to #clusters
%        belongToExpert = vec2ind(gateOut); 
    case 2
       % each data sample has a chance to be classified based on it's
       % probability. for expample, given 2 experts and a datasample with
       % p1=0.6 chance to belong to the 1st expert and p2=0.4 chance to
       % belong to the 2nd expert. so this method will classify this point
       % according to the probability.
       for k=1:size(gateOut,2)
            acumulativeProb = tril(ones(expertCount,expertCount))*gateOut(:,k);
            belongToExpert(1,k) = find(acumulativeProb > rand(1),1);
       end
    case 3
        % the MoE network output is the sum of each expert output times the
        % probability of this expert (this method was taken from
        % jacobs1990 paper).
        for j=1:expertCount
            tempNet = expertsNN{1,j};
            expertsOut(j,:) = tempNet(NNinputs);
        end
        for j=1:size(NNtargets,2)
            netOut(1,j) = (expertsOut(:,j))'*gateOut(:,j);
        end
        targ = NNtargets;
    otherwise
        error('wrong "competetiveFlag", try again');
end

if competetiveFlag==1 || competetiveFlag==2
    if size(NNinputs,2) > 1 % many samples, performance analisys
        % sending the samples to their experts
        for j=1:expertCount
            cluster_i_train_ind{1,j} = find(belongToExpert == j);
            tempNet = expertsNN{1,j};
            out_temp = tempNet(NNinputs(:,cluster_i_train_ind{1,j}));
            targ_temp = NNtargets(:,cluster_i_train_ind{1,j});
            targ = [targ,targ_temp];
            netOut = [netOut,out_temp];
        end
    else
        tempNet = expertsNN{1,belongToExpert};
        netOut = tempNet(NNinputs(:,1));
    end
end

end

