function [ samplNewSymmetric,TargetNewSymmetric ] = matsuoka_symmetric( dataPointsNum,NNinput,NNtarg )
% this function takes a 4-neuron Matsuoka CPG and arranges all of the genes
% in all of the 24 possible CPG structures (how many ways we can 4 neurons
% in a network)

% Note: this function takes all of the matsuoka parameters (without C_i!!!)
% the target can be either the period or the frequency

% inputs:
% 1) 'dataPointsNum' - the amount of points the we want to permutate (total
%                       number of samples is 24*dataPointsNum)
% 2) 'NNinput' - NN inputs
% 3) 'NNtarg' - NN targets

% outputs:
% 1) 'samplNewSymmetric' - the inputs in 24 permutations
% 2) 'TargetNewSymmetric' - the targets for the 24 permutations (repeating
%                           itself!)


% prepare the data points that we want
dataInd4Train = randsample(length(NNtarg),dataPointsNum);
sampl4train = NNinput(:,dataInd4Train);
targ4train = NNtarg(1,dataInd4Train);

mixingVectos = perms([1,2,3,4]); % how many ways we can arrange the Matsuoka network

samplNewSymmetric = sampl4train;
TargetNewSymmetric = targ4train;

for j=1:length(mixingVectos)
    
    mixVector = mixingVectos(j,:);
    newSampl = zeros(size(sampl4train,1),dataPointsNum);
    for i=1:dataPointsNum

        newSampl(1:2,i) = sampl4train(1:2,i);
%         newSampl(3:6,i) = sampl4train(mixVector+2,i); % no C

        mat = [0             ,sampl4train(3,i)   ,sampl4train(4,i) ,sampl4train(5,i);
             sampl4train(6,i),0                  ,sampl4train(7,i) ,sampl4train(8,i);
             sampl4train(9,i),sampl4train(10,i)  ,0                ,sampl4train(11,i);
             sampl4train(12,i),sampl4train(13,i) ,sampl4train(14,i),0               ];

        mat = mat(mixVector,mixVector);

        newSampl(3:14,i) = [mat(1,2);mat(1,3);mat(1,4);...
                            mat(2,1);mat(2,3);mat(2,4);...
                            mat(3,1);mat(3,2);mat(3,4);...
                            mat(4,1);mat(4,2);mat(4,3)];
                        
    end

    samplNewSymmetric = horzcat(samplNewSymmetric,newSampl);
    TargetNewSymmetric = horzcat(TargetNewSymmetric,targ4train);
end

end

