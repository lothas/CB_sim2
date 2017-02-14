function [obj] = paper_MoE_train(obj)
% based on code from the internet: https://goker.wordpress.com/2011/07/01/mixture-of-experts/

% Trains a mixture of experts logistic discriminator with given number of
% experts (all experts are linear NN with no hidden layers)

% TODO: add option to insert initial weights for fast retraining of exsiting networks  

% in this function the data is taken as the transpose. meaning, rows are
% samples and columns are features.
sampl_train = (obj.sampl_train)';
targ_train = (obj.targ_train)';

% the validation samples were not transposed becasuse they are transposed in the
% test function.
sampl_valid = obj.sampl_valid;
targ_valid = obj.targ_valid;

% the test samples were not transposed becasuse they are transposed in the
% test function.
sampl_test = obj.sampl_test;
targ_test = obj.targ_test;

expertCount = obj.expertCount;
learningRate = obj.paper_MoE_out.learningRate;
decay = obj.paper_MoE_out.decay;

inputCount = size(sampl_train,2);
sampleCount = size(sampl_train,1);
dim = inputCount + 1; % hidden neurons num?
outputCount = size(targ_train, 2);
 
% add bias unit to x
sampl_train = [ones(sampleCount,1) sampl_train];
 
% initilize random weights:
[ExpertsWeights,gateWeights] = obj.paper_MoE_init(inputCount,outputCount,expertCount,dim);
 
iters = 1;
errs = zeros(obj.numOfIteretions, 1);
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
    err = obj.paper_MoE_test(sampl_valid, targ_valid, ExpertsWeights, gateWeights,0);
    
    if obj.disp_information
        disp(['at itre num ',num2str(iters),': Error = ', num2str(err),' learning rate is: ',num2str(learningRate)]);
    end
    
    errs(iters, 1) = err;
    % take average error of last five runs
    si = (((iters - 5) >= 0) * (iters - 5)) + 1;
    li = si + 4;
    iters = iters + 1;    
    % check stop condition
    if iters > obj.numOfIteretions
        if obj.disp_information
            fprintf('Max Iterations Reached\n');
        end
        break;
    end
%     if err < 1e-6 || sum( errs(si:li,:) ) < sum( errs(si+1:li+1,:) )
%         fprintf('Error reached minimum\n');
%         break;
%     end
end

[~, netOut_train, gateOut_train] = obj.paper_MoE_test(obj.sampl_train, obj.targ_train, ExpertsWeights, gateWeights,0);
[~, netOut_valid, gateOut_valid] = obj.paper_MoE_test(obj.sampl_valid, obj.targ_valid, ExpertsWeights, gateWeights,0);
[~, netOut_test, gateOut_test] = obj.paper_MoE_test(sampl_test, targ_test, ExpertsWeights, gateWeights,0);

obj.paper_MoE_out.out_from_train = netOut_train;
obj.paper_MoE_out.out_from_valid = netOut_valid;
obj.paper_MoE_out.out_from_test = netOut_test;
obj.paper_MoE_out.gateOut_from_train = gateOut_train;
obj.paper_MoE_out.gateOut_from_valid = gateOut_valid;
obj.paper_MoE_out.gateOut_from_test = gateOut_test;

obj.paper_MoE_out.Moe_perf_over_iter = errs';
obj.paper_MoE_out.ExpertsWeights = ExpertsWeights; 
obj.paper_MoE_out.gateWeights = gateWeights;

end

