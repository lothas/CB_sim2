function [netOut,gateOut,belongToExpert,cluster_i_train_ind] =...
    MoE_testNet(obj,NNinputs,expertsNN,...
    gateNet,MoE_method)
%this function train a "Mixture of Experts" newural networks

%%%%%%%%% TODO: make this code to work with multiple NN outputs! %%%%%%%

% inputs:
% 1) NNinputs - the inputs for the neural networks,
%               matrix of (#inputs x #samples) dimention
% 2) NNtargets - the tarhets for the neural networks,
%                matrix of (#targets x #samples) dimention
% 3) expertsNN - cell array containing the experts and their performance.
% 4) gateNet - the Gate NN.
% 5) MoE_method - *) 'hardCompetetetive' = "winner takes all"
%                 *) 'softCompetetetive' = "chance for everybody"
%                 *) 'collaboration' = (out = expertsOut*gateOut)

% outputs:
% 1) netOut - output from MoE
% 2) gateOut - the gate output, represents the probability of a sample to
%               belong to each expert
% 3) belongToExpert - which expert's cluster the sample belongs.
% 4) cluster_i_train_ind - in case of clustering. this is the indecies
%   specify which samples went to which expert.

expertCount = size(expertsNN,2);
num_of_samples = size(NNinputs,2);

belongToExpert = zeros(1,num_of_samples);
expertsOut = zeros(expertCount,num_of_samples); % Experts output matrix 
cluster_i_train_ind = cell(1,expertCount);

% sending samples (inputs) to gate net:
gateOut = gateNet(NNinputs);

switch MoE_method
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
    case 'collaboration'
        % the MoE network output is the sum of each expert output times the
        % probability of this expert (this method was taken from
        % jacobs1990 paper).
        netOut = zeros(1,num_of_samples);
        for j=1:expertCount
            tempNet = expertsNN{1,j};
            expertsOut(j,:) = tempNet(NNinputs);
        end
        for j=1:num_of_samples
            netOut(1,j) = (expertsOut(:,j))'*gateOut(:,j);
        end
    otherwise
        error('wrong "competetiveFlag", try again');
end

switch MoE_method
    case {1,2}
        if num_of_samples > 1 % many samples, performance analisys
            % sending the samples to their experts
            netOut = zeros(1,num_of_samples);
            for j=1:expertCount
                % check wich sample belongs to which cluster (expert):
                cluster_i_train_ind{1,j} = find(belongToExpert == j);

                % define each expert as a temporary 'net' object:
                tempNet = expertsNN{1,j};

                % make expert outputs and rearrange it according to the
                % original smaples order:
                netOut(1,cluster_i_train_ind{1,j}) = tempNet(NNinputs(:,cluster_i_train_ind{1,j}));
            end
        else
            % if only one sampl was sent for testing:
            tempNet = expertsNN{1,belongToExpert};
            netOut = tempNet(NNinputs(:,1));
        end
end

% if graphGO
%     obj.my_MoE_plot_test_perf(expertCount,netOut,NNtargets,...
%         cluster_i_train_ind,gateOut,competetiveFlag)
% end

end

