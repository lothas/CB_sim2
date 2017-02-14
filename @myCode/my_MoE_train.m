function [obj] = my_MoE_train(obj)
%this function train a "Mixture of Experts" newural networks

sampl_train = obj.sampl_train;
targ_train = obj.targ_train;
sampl_valid = obj.sampl_valid;
targ_valid = obj.targ_valid;
sampl_test = obj.sampl_test;
targ_test = obj.targ_test;

disp('initilizing new MoE:');
tic
obj = obj.My_MoE_init();
disp(['time for init: ',num2str(toc)]);

% defining some constants:
num_of_train_samples = size(sampl_train,2);
numOfIteretions = obj.numOfIteretions;
expertCount = obj.expertCount; % number of "Experts", each Expert is a NN
competetiveFlag = obj.my_MoE_out.competetiveFlag; % our preffered mathod of training:
%                                           if '1'- "winner takes all"
%                                              '2'- "chance for everybody"
%                                              '3'- out = expertsOut * gateOut

gateNet = obj.my_MoE_out.gateNet;
expertsNN = obj.my_MoE_out.expertsNN;

% data storage:
expert_i_GroupSize = zeros(expertCount,numOfIteretions); % num of samples in each expert's cluster
outMat = zeros(expertCount,num_of_train_samples); % Experts output matrix 
errMat = zeros(size(outMat)); % error matrix (each row - targets)
gateNN_perf_vec = zeros(1,numOfIteretions); % gate performance over iteration num
Experts_perf_mat = zeros(expertCount,numOfIteretions); % experts best perf (MSE) over iteration num
gateNet_targ = zeros(expertCount,num_of_train_samples); % target for gate network (note: only one '1' in each coloumn)
emptyGroupIndecator = false(expertCount,numOfIteretions); % count how many time we have an empty cluster
Moe_perf_over_iter = zeros(1,numOfIteretions); % the performance of the entire MoE over #iteretion


disp('start training...');
tic
% Train the experts and the gate NN:

for i=1:numOfIteretions
    disp(['at iteration num: #',num2str(i)]);
    
    % test network to check  validation_error:
    [MoE_out_valid,~,~,~] = obj.my_MoE_testNet(sampl_valid,targ_valid,expertsNN,...
    gateNet,competetiveFlag,0);

    % calc MoE performance to check wether to stop the training:
    [Moe_perf_over_iter(1,i),~] = obj.NN_perf_calc(targ_valid,MoE_out_valid,0,0);
    if Moe_perf_over_iter(1,i) < 0.000000000001 % stopping condition on error
        disp('reached below the desired error');
        break;
    end
    if (i > 11) && (mean(Moe_perf_over_iter(1,(i-10):i)) < 0.00000000001) % stopping condition on error gradient
        disp('reached below the desired error gradient');
        break;
    end
    
    % devide to clusters with the gate network and get the gate output:
    [netOut_train,gateOut,~,cluster_i__ind] = obj.my_MoE_testNet(sampl_train,targ_train,expertsNN,...
    gateNet,competetiveFlag,0);

    % run each expert on the entire data:
    for j=1:expertCount
        tempNet = expertsNN{1,j};
        outMat(j,:) = tempNet(sampl_train);
        errMat(j,:) = outMat(j,:) - targ_train;
    end
    seMat = errMat.^2; % squar error
            
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
                        sampl_train(:,cluster_i__ind{1,j}), targ_train(:,cluster_i__ind{1,j}));
                    % TODO: think what to do if the cluster is empty
                else
                    % if the expert's cluster is empty, train the expert on the 'n'
                    % points with the best probability
                    [~,tempGroup] = sort(gateOut(j,:));
                    best_tempGroup = tempGroup(1,1:1000);
                    [expertsNN{1,j}, expertsNN{2,j}] = train(tempNet,...
                        sampl_train(:,best_tempGroup), targ_train(:,best_tempGroup));
                    emptyGroupIndecator(j,i) = true;
                end
                trExpertPerf_temp = expertsNN{2,j};
                Experts_perf_mat(j,i) = trExpertPerf_temp.best_perf;
            end
            
            [~,best_expert_ind] = min(seMat,[],1);

            % retrain gate network:
            for j=1:expertCount % find targets based on experts perf
                gateNet_targ(j,:) = (best_expert_ind == j);
            end
            [gateNet,gateNet_perf] = train(gateNet,sampl_train,gateNet_targ);
            gateNN_perf_vec(1,i) = gateNet_perf.best_perf;
            
        case 3
            % change the importance of each Target based on the probability
            % of this sample to belong to a certain expert.
            % note: exprimental code, use with caution!
            fh = zeros(expertCount,num_of_train_samples);
            g = gateOut;

            yStar_yi = seMat;
            
            for k=1:size(targ_train,2)
                fh(:,k) = g(:,k) .* exp(-0.5 .* yStar_yi(:,k) );
                fh(:,k) = fh(:,k) ./ sum(fh(:,k),1);
            end
            for j=1:expertCount
                tempNet = expertsNN{1,j};
                errorWeights = {fh(j,:)};                
                [expertsNN{1,j}, expertsNN{2,j}] = train(tempNet,...
                        sampl_train, targ_train,[],[],errorWeights);
            end
%             [fhMax,fhMaxIndex] = max(fh,[],1);
%             % 'fh' becomes the target of the gate network by taking the max
%             % of 'fh' => '1' all other in row => '0'
%             gateNet_targ = zeros(expertCount,num_of_train_samples);
%             gateNet_targ(:,fhMaxIndex) = 1;
%             
%             % 'fh' is also the "errorWeigths" for the training of the gate
%             % network.
%             
%             [gateNet,gateNet_perf] = train(gateNet,sampl_train,gateNet_targ,[],[],fhMax);
%             gateNN_perf_vec(1,i) = gateNet_perf.best_perf;
            
            [gateNet,gateNet_perf] = train(gateNet,sampl_train,fh);
            gateNN_perf_vec(1,i) = gateNet_perf.best_perf;
            
        otherwise
            error('wrong "competetiveFlag", try again');
    end
    
end

   
obj.my_MoE_out.expertsNN = expertsNN;
obj.my_MoE_out.gateNet = gateNet;
obj.my_MoE_out.Moe_perf_over_iter = Moe_perf_over_iter;
obj.my_MoE_out.gateTraniData.gateNN_perf_vec = gateNN_perf_vec;
obj.my_MoE_out.out_from_train = netOut_train;
obj.my_MoE_out.out_from_valid = MoE_out_valid;

% test network to check test error:
    [MoE_out_test,~,~,~] = obj.my_MoE_testNet(sampl_test,targ_test,expertsNN,...
    gateNet,competetiveFlag,0);
obj.my_MoE_out.out_from_test = MoE_out_test;

switch competetiveFlag
    case {1,2}
        obj.my_MoE_out.expertsTrainData.expert_i_GroupSize = expert_i_GroupSize;
        obj.my_MoE_out.expertsTrainData.Experts_perf_mat = Experts_perf_mat;
        obj.my_MoE_out.expertsTrainData.emptyGroupIndecator = emptyGroupIndecator;
    case 3
        obj.my_MoE_out.expertsTrainData.expert_i_GroupSize = [];
        obj.my_MoE_out.expertsTrainData.Experts_perf_mat = [];
        obj.my_MoE_out.expertsTrainData.emptyGroupIndecator = [];
end


disp(['runTime of training: ',num2str(toc)]);

disp([' MoE perf (MSE) on test group: ',num2str(Moe_perf_over_iter(1,end))]);

end

