function [ expertsNN,gateNet,gateNet_perf,...
    expert_i_GroupSize,gateNN_perf_vec,Experts_perf_mat,...
    R_squar,errMSE,emptyGroupIndecator ] = ...
    my_MoE_train(NNinputs,NNtargets,expertCount,...
    numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
    GateHidLayer,GateHidNueron,GraphicsFlag,competetiveFlag)
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
% 8) competetiveFlag - if 'competetive'- than the clustering is being done
%                       by highest P takes all. if not competetive than
%                       each sample as a chance to go to each cluster based
%                       on the probability from the gate network

% outputs:
% 1) expertsNN - cell array containing the experts and their performance. 
% 2) gateNet - the Gate NN.
% 3) gateNet_perf - the training performace of the gateNN
% 4) expert_i_GroupSize - the expert cluster size over interation num
% 5) gateNN_perf_vec - gateNet performance over iteration
% 6) Experts_perf_mat - experts best perf (MSE) over iteration num
% 7) emptyGroupIndecator - bolean matrix, '1'- every iteration that we have
%                           an empty group. for comparing the perf matrix
%                           sudden changes.
% 8) R_squar - R^2 after rechecking the training data
% 9) errMSE - MSE after rechecking the training data

disp(['preparing ' ,num2str(expertCount),' Experts with ',...
    num2str(ExpertHidLayer),' hidden layers and ',num2str(GateHidNueron),...
    ' neurons in each layer...']);
disp(['preparing gate NN with ',...
    num2str(GateHidLayer),' hidden layers and ',num2str(GateHidNueron),...
    ' neurons in each layer...']);

expert_i_GroupSize = zeros(expertCount,numOfIteretions); % num of samples in each expert's cluster
cluster_i__ind = cell(1,expertCount); % indecies vector of samples for each expert
outMat = zeros(expertCount,size(NNinputs,2)); % Experts output matrix 
errMat = zeros(size(outMat)); % error matrix (each row - targets)
gateNN_perf_vec = zeros(1,numOfIteretions); % gate performance over iteration num
Experts_perf_mat = zeros(expertCount,numOfIteretions); % experts best perf (MSE) over iteration num
gateNet_targ = zeros(expertCount,size(NNinputs,2)); % target for gate network (note: only one '1' in each coloumn)
emptyGroupIndecator = false(expertCount,numOfIteretions); % count how many time we have an empty cluster
expertsNN = cell(2,expertCount);


% define experts:
% define each expert as a NN with 'ExpertHidLayer' with 'ExpertHidNueron'
% hidden neurons each.
expertsBuildVector = ExpertHidNueron*ones(1,ExpertHidLayer);
for j=1:expertCount
    expertsNN{1,j} = feedforwardnet(expertsBuildVector);
    expertsNN{1,j}.trainParam.showWindow = 0; % dont show training window
    expertsNN{1,j}.trainParam.epochs = maxEphocs;
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
gateNet.trainParam.showWindow = 0;
gateNet = train(gateNet,NNinputs,gateNet_targ);

disp('start training...');
%% Train the experts and the gate NN:
tic
for i=1:numOfIteretions
    disp(['at iteration num: #',num2str(i)]);
    
    % devide to clusters with the gate network:
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
    % clustering to different experts
    for j=1:expertCount
        cluster_i__ind{1,j} = find(gateOut_inx == j);
        expert_i_GroupSize(j,i) = length(cluster_i__ind{1,j}); % check the size of each cluster
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
end

% check the final clustering after training (used for MoE train perf checking)
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

disp('runTime of training:');
toc

if GraphicsFlag % some performance graphs
    for j=1:expertCount
       figure;
       groupInd = cluster_i__ind{1,j};
       tempNet = expertsNN{1,j};
       Outputs = tempNet(NNinputs);
       trOut = Outputs(:,groupInd);
       trTarg = NNtargets(groupInd);
       plotregression(trTarg,trOut,'Train');
       clear groupInd tempNet Outputs trOut
    end

    expertsNames = cell(1,expertCount);
    for j=1:expertCount
        expertsNames{1,j} = ['#',num2str(j),' expert'];
    end
    figure;
    plot(1:numOfIteretions,expert_i_GroupSize); hold on;
    xlabel('#iteretion');   ylabel('group size [#points]');
    legend(expertsNames);
    
    bestExpertsInx = [gateOut_inx;-best_expert_ind];
    disp(['the number of train points correctly classify: ',...
        num2str(length(find(~(sum(bestExpertsInx,1))))),...
        ' out of ',num2str(length(bestExpertsInx))]);
end

%% check train results:
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

% calc R^2
err = targ_train-outM_train;
errVar = var(err,0,2);
inputVar = var(targ_train,0,2);

R_squar = 1-(errVar/inputVar);
errMSE = immse(outM_train,targ_train);

if GraphicsFlag % some performance graphs
    figure;
    plotregression(targ_train,outM_train,'Train');
    disp('checking the training performance:');
    disp(['the R^2 is: ',num2str(R_squar)]);
    disp(['the MSE is: ',num2str(errMSE)]);
end
end

