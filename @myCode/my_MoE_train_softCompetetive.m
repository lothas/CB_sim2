function [obj] = my_MoE_train_softCompetetive(obj)
%this function train a "Mixture of Experts" newural networks
% % NOTE: same for the "my_MoE_train.m" but only for competetiveFlag=2 and
% with better visualization.

switch obj.expertCount
    case {2,3} % in case of small number of expert, make colors clear:
        colors = [1,0,0;0,1,0;0,0,1];
    otherwise
        colors = rand(obj.expertCount,3);
end

fig1 = figure;
ax1 = subplot(2,2,1);% axes('Position',[0.1 0.1 0.7 0.7]);
ax2 = subplot(2,2,2);
ax3 = subplot(2,2,3);
ax4 = subplot(2,2,4);

ax1.NextPlot= 'add';
ax1.XLabel.String = '#iteration';
ax1.YLabel.String = 'MSE of MoE net';   
ax1.Title.String = 'MSE perf over iteration';
ax1.XGrid = 'on';   ax1.XMinorGrid  = 'on';
ax1.YGrid = 'on';   ax1.YMinorGrid  = 'on';

ax2.NextPlot= 'add';
ax2.XLabel.String = '#iteration';
ax2.YLabel.String = 'gate performance MSE of g-f_{h}';
ax2.Title.String = 'gate perf over iteration';
ax2.XGrid = 'on';   ax2.XMinorGrid  = 'on';
ax2.YGrid = 'on';   ax2.YMinorGrid  = 'on';

ax3.NextPlot= 'add';
ax3.XLabel.String='target'; ax3.YLabel.String='output';
ax3.Title.String={'regression graph with different color for every dominant expert';'empty circle mean g<0.5'};
line([0,1],[0,1],'color',[0,0,0],'linewidth',2,'Parent',ax3);

ax4.NextPlot= 'add';
ax4.XLabel.String='sample Num'; 
ax4.YLabel.String='gate output [Prob]';
ax4.Title.String={'The probability of each sample to belong to each expert',...
    'each expert is a different color'};

sampl_train = obj.sampl_train;
targ_train = obj.targ_train;
sampl_valid = obj.sampl_valid;
targ_valid = obj.targ_valid;
sampl_test = obj.sampl_test;
targ_test = obj.targ_test;

% take only few points for the regression plot:
numToShow = 100;
if size(targ_test,2) > numToShow
    randSampl_ind = randsample(size(targ_test,2),numToShow);
    randSampl_ind = randSampl_ind';
else
    numToShow = size(targ_test,2);
    randSampl_ind = 1:size(targ_test,2);
end

% Initializing the Experts:
if obj.disp_information
    disp('initilizing new MoE:');
    tic
end
obj = obj.My_MoE_init();
if obj.disp_information
    disp(['time for init: ',num2str(toc)]);
end

% defining some constants:
num_of_train_samples = size(sampl_train,2);
numOfIteretions = obj.numOfIteretions;
expertCount = obj.expertCount; % number of "Experts", each Expert is a NN
competetiveFlag = 2;

gateNet = obj.my_MoE_out.gateNet;
expertsNN = obj.my_MoE_out.expertsNN;

% data storage:
outMat = zeros(expertCount,num_of_train_samples); % Experts output matrix 
errMat = zeros(size(outMat)); % error matrix (each row - targets)
gateNN_perf_vec = zeros(1,numOfIteretions); % gate performance over iteration num
Moe_perf_over_iter = zeros(1,numOfIteretions); % the performance of the entire MoE over #iteretion
Experts_perf_mat = zeros(expertCount,numOfIteretions); % experts best perf (MSE) over iteration num
emptyGroupIndecator = false(expertCount,numOfIteretions); % count how many time we have an empty cluster

if obj.disp_information
    disp('start training...');
    tic
end

% Train the experts and the gate NN:

for i=1:numOfIteretions
            
    % update the graphs:
    drawnow
    
    if obj.disp_information
        disp(['at iteration num: #',num2str(i)]);
    end
    
    % test network (on various group) to check  validation_error:
    [netOut_train,gateOut,~,cluster_i__ind] = obj.my_MoE_testNet(sampl_train,targ_train,expertsNN,...
    gateNet,competetiveFlag,0);
    [MoE_out_valid,~,~,~] = obj.my_MoE_testNet(sampl_valid,targ_valid,expertsNN,...
    gateNet,competetiveFlag,0);
    [MoE_out_test,gateOut_test,~,~] = obj.my_MoE_testNet(sampl_test,targ_test,expertsNN,...
    gateNet,competetiveFlag,0);

    [MSE_train,~] = obj.NN_perf_calc(targ_train,netOut_train,0,0);
    [MSE_valid,~] = obj.NN_perf_calc(targ_valid,MoE_out_valid,0,0);
    [MSE_test,~] = obj.NN_perf_calc(targ_test,MoE_out_test,0,0);
    
    % plot the MoE performance
    plot(ax1,i,MSE_train,'b-o');
    plot(ax1,i,MSE_valid,'g-o');
    plot(ax1,i,MSE_test,'r-o');
    
    % plot and update the regression plot:
    [g_max,g_max_ind] = max(gateOut_test,[],1);
    ouputs = MoE_out_test;
    targets = targ_test;
    cla(ax3);
    for k=1:obj.expertCount
        for n=1:numToShow
            ind = randSampl_ind(1,n);
            if g_max_ind(1,ind) == k
                if g_max(1,ind) > 0.5
                    plot(ax3,targets(1,ind),ouputs(1,ind),'k-o','MarkerFaceColor',colors(k,:));
                else
                    plot(ax3,targets(1,ind),ouputs(1,ind),'k-o');
                end
            end
        end
    end
    
    cla(ax4);
    bar(ax4,(gateOut_test(:,randSampl_ind))','stacked');
    
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
    
   % check the cluster size of each expert: 
   for j=1:expertCount % check the size of each cluster
        expert_i_GroupSize(j,i) = length(cluster_i__ind{1,j}); 
   end
   
    % run each expert on the entire data:
    for j=1:expertCount
        tempNet = expertsNN{1,j};
        outMat(j,:) = tempNet(sampl_train);
        errMat(j,:) = outMat(j,:) - targ_train;
    end
    seMat = errMat.^2; % squar error
    % initilize f_h
    fh = zeros(expertCount,num_of_train_samples);
    g = gateOut;
    yStar_yi = seMat;
    % calc f_h from the equation in Jacobs1990 paper
    for k=1:size(targ_train,2)
        fh(:,k) = g(:,k) .* exp(-0.5 .* yStar_yi(:,k) );
        fh(:,k) = fh(:,k) ./ sum(fh(:,k),1);
    end
   
    % training the experts:
    %   train each expert only on the samples which belongs to it.
    %      also train the experts with weights given to samples by f_h 
    for j=1:expertCount
        tempNet = expertsNN{1,j};
        errorWeights_all0 = {fh(j,:)};
        errorWeights_all = errorWeights_all0{1,1};
        if expert_i_GroupSize(j,i) > 0 % only train if the cluster is not empty
            errorWeights = errorWeights_all(1,cluster_i__ind{1,j});
            [expertsNN{1,j}, expertsNN{2,j}] = train(tempNet,...
                sampl_train(:,cluster_i__ind{1,j}),...
                targ_train(:,cluster_i__ind{1,j}),...
                [],[],errorWeights);
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
            

    % train the gate using f_h as targets:
    %       (minimize MSE between 'g' and 'f_h')
    [gateNet,gateNet_perf] = train(gateNet,sampl_train,fh);
    plot(ax2,i,gateNet_perf.best_perf,'k-o');
    gateNN_perf_vec(1,i) = gateNet_perf.best_perf;
    
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
[Moe_MSE_on_test_temp,Rsquar] = obj.NN_perf_calc(targ_test,MoE_out_test,0,0);


obj.my_MoE_out.expertsTrainData.expert_i_GroupSize = expert_i_GroupSize;
obj.my_MoE_out.expertsTrainData.Experts_perf_mat = Experts_perf_mat;
obj.my_MoE_out.expertsTrainData.emptyGroupIndecator = emptyGroupIndecator;

obj.my_MoE_out.Moe_MSE_on_test = Moe_MSE_on_test_temp;
obj.my_MoE_out.RsquarTest = Rsquar;
if obj.disp_information
    disp(['runTime of training: ',num2str(toc)]);
    disp([' MoE perf (MSE) on test group: ',num2str(Moe_MSE_on_test_temp)]);
end

end

