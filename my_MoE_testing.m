function [R_squar,errMSE] = my_MoE_testing(NNinputs,NNtargets,expertsNN,...
    gateNet,GraphicsFlag,competetiveFlag)
%this function train a "Mixture of Experts" newural networks

% inputs:
% 1) NNinputs - the inputs for the neural networks,
%               matrix of (#inputs x #samples) dimention
% 2) NNtargets - the tarhets for the neural networks,
%                matrix of (#targets x #samples) dimention
% 3) expertsNN - cell array containing the experts and their performance.
% 4) gateNet - the Gate NN.
% 5) competetiveFlag - if 'competetive'- than the clustering is being done
%                       by highest P takes all. if not competetive than
%                       each sample as a chance to go to each cluster based
%                       on the probability from the gate network

% outputs:
% 1) R_squar - R^2 after rechecking the training data
% 2) errMSE - MSE after rechecking the training data

expertCount = size(expertsNN,2);

disp('start testing...');
tic

% check the final clustering after training (used for MoE train perf checking)
gateOut_inx = zeros(1,size(NNinputs,2));
gateOut = gateNet(NNinputs);
if competetiveFlag
    gateOut_inx = vec2ind(gateOut); % indecies of points to #clusters
else
   % cluster by probability and not by highest chance
   for k=1:size(gateOut,2)
        acumulativeProb = tril(ones(expertCount,expertCount))*gateOut(:,k);
        gateOut_inx(1,k) = find(acumulativeProb > rand(1),1);
   end
end
disp('total runTime of testing:');
toc

% run each expert on the entire data: check how many were clusify correctly
% based on minimum squar error.
outMat = zeros(expertCount,size(NNinputs,2)); % Experts output matrix 
errMat = zeros(size(outMat)); % error matrix (each row - targets)
for j=1:expertCount
    tempNet = expertsNN{1,j};
    outMat(j,:) = tempNet(NNinputs);
    errMat(j,:) = outMat(j,:) - NNtargets;
end
seMat = errMat.^2; % squar error
[~,best_expert_ind] = min(seMat,[],1); 
bestExpertsInx = [gateOut_inx;-best_expert_ind];
disp(['the number of train points correctly classify: ',...
    num2str(length(find(~(sum(bestExpertsInx,1))))),...
    ' out of ',num2str(length(bestExpertsInx))]);

%% check test results:
cluster_i_train_ind = cell(1,expertCount);
targ_train = [];
outM_train = [];
for j=1:expertCount
    cluster_i_train_ind{1,j} = find(gateOut_inx == j);
    tempNet = expertsNN{1,j};
    out_train_temp = tempNet(NNinputs(:,cluster_i_train_ind{1,j}));
    targ_temp = NNtargets(:,cluster_i_train_ind{1,j});
    targ_train = [targ_train,targ_temp];
    outM_train = [outM_train,out_train_temp];
    clear targ_temp out_train_temp
end

% show performance:
[errMSE,R_squar] = NN_perf_calc(Targets,NNoutput,1,GraphicsFlag );
end

