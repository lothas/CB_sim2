function [ExpertsWeights, gateWeights,errs,trainingIds,testingIds] = paper_MoE_train(sampl, targ, expertCount, maxIter, learningRate, decay)
% based on code from the internet: https://goker.wordpress.com/2011/07/01/mixture-of-experts/

% Trains a mixture of experts logistic discriminator with given number of
% experts (all experts are linear NN)

% sampl: training data, samples in rows
% targ: output values in each row for regression problems
% expertCount: number of experts (logistic discriminators) to use
% maxIter: a stopping cond for the loop
% learningRate: learning rate for NN training
% decay: decay of the learning rate 

% TODO: add option to insert initial weights for fast retraining of exsiting networks  

[sampl_train,targ_train,trainingIds,...
    sampl_test,targ_test,testingIds ] = divide_train_and_test(sampl,targ,0.8);

sampl_train = sampl_train';
targ_train = targ_train'; 

sampleCount = size(sampl_train,1);
dim = size(sampl_train,2) + 1; % hidden neurons num?
outputCount = size(targ_train, 2);
 
% add bias unit to x
sampl_train = [ones(sampleCount,1) sampl_train];
 
% initialize parameters:
% initilize random weights by taking a random number from -1 to 1 devide by number of weights.
% gateNet weights-
gateWeights = (2*rand(expertCount, dim)-1)/size(sampl_train,2); 
% experts weights-
ExpertsWeights = (2*rand(expertCount,dim, outputCount)-1)/size(sampl_train,2);
 
iters = 1;
errs = zeros(maxIter, 1);
while true
    % choose next training instance randomly
    trSeq = randperm(sampleCount);
    for i=1:sampleCount
        k = trSeq(i);
        train_sampl = sampl_train(k,:);
        train_targ = targ_train(k,:);
        % calculate intermediate values
        g = (exp(gateWeights*train_sampl'))';
        g = g ./ sum(g);
        ExpertOut = zeros(outputCount, expertCount); % ExpertOut-experts out

        for j=1:outputCount % calculate output from the experts
            ExpertOut(j,:,:) = (ExpertsWeights(:,:,j)*train_sampl')';  
        end        
        
        % output for each output dimension
%         yi = (ExpertOut*g')';
        
        % update parameters
        fh = exp(sum( (repmat( train_targ', 1, expertCount )  - ExpertOut) .^2, 1 ) .* -0.5) .* g;
    
        fh = fh ./ sum(fh);
        
        % calculate delta ExpertsWeights and gateWeights for each output unit
        for r=1:expertCount
            for j=1:outputCount
                dv = learningRate .* ( train_targ(1,j) - ExpertOut(j, r)) * fh(1, r) * train_sampl;
                ExpertsWeights(r, :, j) = ExpertsWeights(r, :, j) + dv;
            end

            % calculate delta gateWeights
            dm = learningRate .* ( fh(1, r) - g(1, r) ) * train_sampl;
            gateWeights( r, : ) = gateWeights( r, : ) + dm;
        end        
    end
    learningRate = learningRate * decay;
    
    % calculate training set error
    err = paper_MoE_test(sampl_test, targ_test, ExpertsWeights, gateWeights,0);
    
    disp(['at itre num ',num2str(iters),': Error = ', num2str(err),' learning rate is: ',num2str(learningRate)]);
    errs(iters, 1) = err;
    % take average error of last five runs
    si = (((iters - 5) >= 0) * (iters - 5)) + 1;
    li = si + 4;
    iters = iters + 1;    
    % check stop condition
    if iters > maxIter 
        fprintf('Max Iterations Reached\n');
        break;
    end
%     if err < 1e-6 || sum( errs(si:li,:) ) < sum( errs(si+1:li+1,:) )
%         fprintf('Error reached minimum\n');
%         break;
%     end
end

errs = errs';

end

