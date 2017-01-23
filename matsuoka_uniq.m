function [ samplNewUniq,TargetNewUniq ] = matsuoka_uniq( dataPointsNum,NNinput,NNtarg )
% this function takes a 4-neuron Matsuoka CPG and arranges all of the genes
% in a unique way based on the largest weight for every matsuoka neuron

% Note: this function takes all of the matsuoka parameters (without C_i!!!)
% the target can be either the period or the frequency

% inputs:
% 1) 'dataPointsNum' - the amount of points the we want to use
% 2) 'NNinput' - NN inputs [tau, b, 12 x w_ij] (no c_i !!)
% 3) 'NNtarg' - NN targets

% outputs:
% 1) 'samplNewUniq' - the inputs in a uniq way
% 2) 'TargetNewUniq' - the targets


% prepare the data points that we want
dataInd4Train = randsample(length(NNtarg),dataPointsNum);
sampl4train = NNinput(:,dataInd4Train);
targ4train = NNtarg(1,dataInd4Train);
newSampl = zeros(size(sampl4train,1),dataPointsNum);

for i=1:dataPointsNum

    newSampl(1:2,i) = sampl4train(1:2,i);

    mat = [0             ,sampl4train(3,i)  ,sampl4train(4,i) ,sampl4train(5,i);
        sampl4train(6,i),0                 ,sampl4train(7,i),sampl4train(8,i);
        sampl4train(9,i),sampl4train(10,i) ,0                ,sampl4train(11,i);
        sampl4train(12,i),sampl4train(13,i) ,sampl4train(14,i),0               ];
    
    biggestWeights = max(mat,[],1);
    [~,sortedWeightsIND] = sort(biggestWeights,2);
%     newSampl(2+(1:4),i) = sampl4train(sortedWeightsIND+2,i); % no C_i
    
    mat = mat(sortedWeightsIND',sortedWeightsIND');

    newSampl(3:14,i) = [mat(1,2);mat(1,3);mat(1,4);...
                        mat(2,1);mat(2,3);mat(2,4);...
                        mat(3,1);mat(3,2);mat(3,4);...
                        mat(4,1);mat(4,2);mat(4,3)];


end
    samplNewUniq = newSampl;
    TargetNewUniq = targ4train;
end

