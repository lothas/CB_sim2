function [ expertsNN,gateNet,expert_i_GroupSize,gateNN_perf_vec,Experts_perf_mat,Moe_perf_over_iter,emptyGroupIndecator ] = ...
    my_MoE_train(NNinputs,NNtargets,expertCount,numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
                GateHidLayer,GateHidNueron,competetiveFlag)
%this function train a "Mixture of Experts" newural networks

% inputs:
% 1) NNinputs - the inputs for the neural networks,
%               matrix of (#inputs x #samples) dimention
% 2) NNtargets - the tarhets for the neural networks,
%                matrix of (#targets x #samples) dimention
% 3) expertCount - number of "Experts" each Expert is a NN.
% 4) numOfIteretions - num of interation of Experts training
% 5) maxEphocs - max num of ephocs in each expert training.
% 6) ExpertHidLayer - num of hidden layers for each expert.
% 7) ExpertHidNueron - num of hidden neurons in each layer.
% 6) GateHidLayer - num of hidden layers for the gate NN.
% 7) GateHidNueron - num of hidden neurons in each layer in gate NN.
% 8) competetiveFlag - if '1'- "winner takes all"
%                         '2'- "chance for everybody"
%                         '3'- out = expertsOut * gateOut

% outputs:
% 1) expertsNN - cell array containing the experts and their performance. 
% 2) gateNet - the Gate NN.
% 3) expert_i_GroupSize - the expert cluster size over interation num
% 4) gateNN_perf_vec - gateNet performance over iteration
% 5) Experts_perf_mat - experts best perf (MSE) over iteration num
% 6) Moe_perf_over_iter - the performance of the entire MoE over #iteretion
% 7) emptyGroupIndecator - bolean matrix, '1'- every iteration that we have
%                           an empty group. for comparing the perf matrix
%                           sudden changes.


%%
disp(['preparing ' ,num2str(expertCount),' Experts with ',...
    num2str(ExpertHidLayer),' hidden layers and ',num2str(GateHidNueron),...
    ' neurons in each layer...']);
disp(['preparing gate NN with ',...
    num2str(GateHidLayer),' hidden layers and ',num2str(GateHidNueron),...
    ' neurons in each layer...']);

expert_i_GroupSize = zeros(expertCount,numOfIteretions); % num of samples in each expert's cluster
outMat = zeros(expertCount,size(NNinputs,2)); % Experts output matrix 
errMat = zeros(size(outMat)); % error matrix (each row - targets)
gateNN_perf_vec = zeros(1,numOfIteretions); % gate performance over iteration num
Experts_perf_mat = zeros(expertCount,numOfIteretions); % experts best perf (MSE) over iteration num
gateNet_targ = zeros(expertCount,size(NNinputs,2)); % target for gate network (note: only one '1' in each coloumn)
emptyGroupIndecator = false(expertCount,numOfIteretions); % count how many time we have an empty cluster
expertsNN = cell(2,expertCount);
Moe_perf_over_iter = zeros(1,numOfIteretions); % the performance of the entire MoE over #iteretion

tic
% define experts:
% define each expert as a NN with 'ExpertHidLayer' with 'ExpertHidNueron'
% hidden neurons each.
expertsBuildVector = ExpertHidNueron*ones(1,ExpertHidLayer);
for j=1:expertCount
    expertsNN{1,j} = feedforwardnet(expertsBuildVector);
    expertsNN{1,j}.trainParam.showWindow = 0; % dont show training window
    expertsNN{1,j}.trainParam.epochs = maxEphocs;
    expertsNN{1,j}.divideMode = 'none'; % all data to training
end

% initial training (to initilazied the NN)
% train each experts on a random sample of 1000 points to initialze the
% weights.
for j=1:expertCount
    initTrain = randsample(1:size(NNinputs,2),1000);
    [expertsNN{1,j}, expertsNN{2,j}] = train(expertsNN{1,j}, NNinputs(:,initTrain), NNtargets(:,initTrain));
end

% define and initialize gate Network:
% define the GateNet as a NN with 'ExpertHidLayer' with 'ExpertHidNueron'
% hidden neurons each.
gateBuildVector = GateHidNueron*ones(1,GateHidLayer);
for j=1:expertCount % run the data throught the experts to get initial clustering
    tempNet = expertsNN{1,j};
    outMat(j,:) = tempNet(NNinputs);
    errMat(j,:) = outMat(j,:) - NNtargets;
end
seMat = errMat.^2; % squar error
[~,best_expert_ind] = min(seMat,[],1);
for j=1:expertCount % find targets based on experts perf
    gateNet_targ(j,:) = (best_expert_ind == j);
end
gateNet = patternnet(gateBuildVector);
gateNet.performParam.regularization = 0; % 0.1;
gateNet.trainParam.showWindow = 0;
gateNet = train(gateNet,NNinputs,gateNet_targ);

disp('initialization time:'); 
toc

disp('start training...');
%% Train the experts and the gate NN:
tic
for i=1:numOfIteretions
    disp(['at iteration num: #',num2str(i)]);
    
    % devide to clusters with the gate network:
    [MoE_out,gateOut,MoE_targ,~,cluster_i__ind] = my_MoE_testNet(NNinputs,NNtargets,expertsNN,...
    gateNet,competetiveFlag);
    
    % calc MoE performance:
    [Moe_perf_over_iter(1,i),~] = NN_perf_calc(MoE_targ,MoE_out,0,0);
    if Moe_perf_over_iter(1,i) < 0.001 % stopping condition on error
        disp('reached below the desired error');
        break;
    end
    if (i > 11) && (mean(Moe_perf_over_iter(1,(i-10):i)) < 0.00001) % stopping condition on error gradient
        disp('reached below the desired error gradient');
        break;
    end
    
    
    switch competetiveFlag
        case {1,2}
           for j=1:expertCount % check the size of each cluster
                expert_i_GroupSize(j,i) = length(cluster_i__ind{1,j}); 
           end
           
            % training the experts;
            for j=1:expertCount
                tempNet = expertsNN{1,j};
                if expert_i_GroupSize(j,i) > 0 % only train if the cluster is not empty
                    [expertsNN{1,j}, expertsNN{2,j}] = train(tempNet,...
                        NNinputs(:,cluster_i__ind{1,j}), NNtargets(:,cluster_i__ind{1,j}));
                    % TODO: think what to do if the cluster is empty
                else
                    % if the expert's cluster is empty, train the expert on the 'n'
                    % points with the best probability
                    [~,tempGroup] = sort(gateOut(j,:));
                    best_tempGroup = tempGroup(1,1:1000);
                    [expertsNN{1,j}, expertsNN{2,j}] = train(tempNet,...
                        NNinputs(:,best_tempGroup), NNtargets(:,best_tempGroup));
                    emptyGroupIndecator(j,i) = true;
                end
                trExpertPerf_temp = expertsNN{2,j};
                Experts_perf_mat(j,i) = trExpertPerf_temp.best_perf;
            end
            % run each expert on the entire data:
            for j=1:expertCount
                tempNet = expertsNN{1,j};
                outMat(j,:) = tempNet(NNinputs);
                errMat(j,:) = outMat(j,:) - NNtargets;
            end
            seMat = errMat.^2; % squar error
            [~,best_expert_ind] = min(seMat,[],1);

            % retrain gate network:
            for j=1:expertCount % find targets based on experts perf
                gateNet_targ(j,:) = (best_expert_ind == j);
            end
            [gateNet,gateNet_perf] = train(gateNet,NNinputs,gateNet_targ);
            gateNN_perf_vec(1,i) = gateNet_perf.best_perf;
            
        case 3
            % change the importance of each Target based on the probability
            % of this sample to belong to a certain expert.
            % note: exprimental code, use with caution!
            ExpertsOuts = zeros(expertCount,size(NNtargets,2));
            fh = zeros(expertCount,size(NNtargets,2));
            g = gateOut;
            for j=1:expertCount
                tempNet = expertsNN{1,j};
                ExpertsOuts(j,:) = tempNet(NNinputs);
            end
            yStar_yi = (repmat(NNtargets,expertCount,1) - ExpertsOuts).^2;
            
            for k=1:size(NNtargets,2)
                fh(:,k) = g(:,k) .* exp(-0.5 .* yStar_yi(:,k) );
                fh(:,k) = fh(:,k) ./ sum(fh(:,k),1);
            end
            for j=1:expertCount
                tempNet = expertsNN{1,j};
                errorWeights = {fh(j,:)};                
                [expertsNN{1,j}, expertsNN{2,j}] = train(tempNet,...
                        NNinputs, NNtargets,[],[],errorWeights);
            end
            [fhMax,fhMaxIndex] = max(fh,[],1);
            % 'fh' becomes the target of the gate network by taking the max
            % of 'fh' => '1' all other in row => '0'
            gateNet_targ = zeros(expertCount,size(NNinputs,2));
            gateNet_targ(:,fhMaxIndex) = 1;
            
            % 'fh' is also the "errorWeigths" for the training of the gate
            % network.
            
            [gateNet,gateNet_perf] = train(gateNet,NNinputs,gateNet_targ,[],[],fhMax);
            gateNN_perf_vec(1,i) = gateNet_perf.best_perf;
            
        otherwise
            error('wrong "competetiveFlag", try again');
    end
    
end

disp('runTime of training:');
toc

end

