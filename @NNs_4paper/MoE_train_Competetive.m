function [obj] = MoE_train_Competetive(obj,type,...
    numOfIteretions,architecture,inputs,targets)
%this function train a "Mixture of Experts" newural networks
% % NOTE: same for the "my_MoE_train.m" but only for competetiveFlag=2 and
% with better visualization.

% type - 1) 'soft' - for "softcompetetive" evry point has a chance to go to
%                   each expert.
%       2) 'hard' - the point goes to the best expert ('highest
%                   probability).

obj.MoE.numOfIteretions = numOfIteretions;
obj.MoE.MoE_method = type;

samplesNum = size(targets,2);

numToShow = 100; % the #samples to show in the online regression plot.

% initilize graphs:
if obj.disp_information
    fig1 = figure;
    set(fig1, 'Position', get(0, 'Screensize')); % set the figure on full screen size

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
    ax4.XLim = [0 numToShow];
    ax4.XLimMode = 'manual';
    ax4.YLim = [0 1];
    ax4.YLimMode = 'manual';
end

[trainInd,valInd,testInd] = dividerand(samplesNum,0.75,0.15,0.15);

inputs_train = inputs(:,trainInd);
targ_train = targets(:,trainInd);
inputs_valid = inputs(:,valInd);
targ_valid = targets(:,valInd);
inputs_test = inputs(:,testInd);
targ_test = targets(:,testInd);

trainSamplesNum = size(targ_train,2);
validSamplesNum = size(targ_valid,2);
testSamplesNum = size(targ_test,2);

% take only few points for the regression plot:
if testSamplesNum > numToShow
    randSampl_ind = randsample(testSamplesNum,numToShow);
    randSampl_ind = randSampl_ind';
else
    numToShow = testSamplesNum;
    randSampl_ind = 1:testSamplesNum;
end

% Initializing the Experts:
MaxEpochs = 5;
expertCount = 3;
obj = obj.MoE_init(expertCount,architecture,MaxEpochs,...
    inputs,targets);

% defining some constants:
expertCount = obj.MoE.expertCount; % number of "Experts", each Expert is a NN
gateNet = obj.MoE.gateNet;
expertsNN = obj.MoE.expertsNN;

% data storage:
gateNN_perf_vec = zeros(1,numOfIteretions); % gate performance over iteration num
Moe_perf_over_iter = zeros(1,numOfIteretions); % the performance of the entire MoE over #iteretion
Experts_perf_mat = zeros(expertCount,numOfIteretions); % experts best perf (MSE) over iteration num
emptyGroupIndecator = false(expertCount,numOfIteretions); % count how many time we have an empty cluster

% some more variables
gateOut_old = zeros(expertCount,testSamplesNum);
g_changes = zeros(numOfIteretions,testSamplesNum);

if obj.disp_information
    disp('start training...');
    tic
end

% Train the experts and the gate NN:

for i=1:numOfIteretions
            
    % update the graphs:
    drawnow
    
    % test network (on various group) to check  validation_error:
    [netOut_train,gateOut,~,cluster_i__ind] = ...
        obj.MoE_testNet(inputs_train,expertsNN,...
        gateNet,type);
    [MoE_out_valid,~,~,~] = ...
        obj.MoE_testNet(inputs_valid,expertsNN,...
        gateNet,type);
    [MoE_out_test,gateOut_test,~,~] = ...
        obj.MoE_testNet(inputs_test,expertsNN,...
        gateNet,type);

    [MSE_train,~] = obj.MoE_perf_calc(targ_train,netOut_train,0,0);
    [MSE_valid,~] = obj.MoE_perf_calc(targ_valid,MoE_out_valid,0,0);
    [MSE_test,~] = obj.MoE_perf_calc(targ_test,MoE_out_test,0,0);
    
    % update some of the graphs
    if obj.disp_information
        % plot the MoE performance
        plot(ax1,i,MSE_train,'b','Marker','o');
        plot(ax1,i,MSE_valid,'g','Marker','o');
        plot(ax1,i,MSE_test,'r','Marker','o');
        ax1.Title.String = ['MSE perf over iteration, perf_{mse} = ',...
            num2str(MSE_test)];
    
    % plot and update the regression plot:
        [g_max,g_max_ind] = max(gateOut_test,[],1);
        ouputs = MoE_out_test;
        cla(ax3);
        for k=1:expertCount
            for n=1:numToShow
                line([0,1],[0,1],'color',[0,0,0],'linewidth',2,'Parent',ax3);
                ind = randSampl_ind(1,n);
                if g_max_ind(1,ind) == k
                    if g_max(1,ind) > 0.5
                        plot(ax3,targ_test(1,ind),ouputs(1,ind),'k-o','MarkerFaceColor',obj.MoE.colors(k,:));
                    else
                        plot(ax3,targ_test(1,ind),ouputs(1,ind),'k-o');
                    end
                end
            end
        end
    
    % bar graph with the gate output of some random points:
        cla(ax4);
        bplot = bar(ax4,(gateOut_test(:,randSampl_ind))','stacked');
        for k=1:expertCount
            set(bplot(k),'facecolor',obj.MoE.colors(k,:))
        end
        AXlegend = legend(bplot, obj.MoE.legendNames,...
            'Location','northwestoutside','FontSize',8);
    end
    
    % calc MoE performance to check wether to stop the training:
    [Moe_perf_over_iter(1,i),~] =...
        obj.MoE_perf_calc(targ_valid,MoE_out_valid,0,0);
    
    perf_Stop_cond = 0.001;
    if Moe_perf_over_iter(1,i) < perf_Stop_cond % stopping condition on error
        disp('reached below the desired error');
        break;
    end
    
    gradient_stop = 0.00001;
    if (i > 11) && (mean(Moe_perf_over_iter(1,(i-10):i)) < gradient_stop) % stopping condition on error gradient
        disp('reached below the desired error gradient');
        break;
        % TODO: use 'diff' to calculate the gradient.
    end
    
   % check the cluster size of each expert: 
   for j=1:expertCount % check the size of each cluster
        expert_i_GroupSize(j,i) = length(cluster_i__ind{1,j}); 
   end
   
   % run each expert on the entire data and calc "fh":
   fh = obj.MoE_calc_fh(inputs_train,targ_train,expertsNN,gateOut); % ( g = gateOut )

    % training the experts:
    %   train each expert only on the samples which belongs to it.
    %      also train the experts with weights given to samples by f_h 
    for j=1:expertCount
        tempNet = expertsNN{1,j};
        errorWeights_all = fh(j,:);
        % only train the experts if they have data points in their group:
        %    NOTE: with soft competetive the expert's cluster can't be empty
        %           (it's not likely)
        if ~isempty(cluster_i__ind{1,j})
            errorWeights = errorWeights_all(1,cluster_i__ind{1,j});
            [tempNet, tempNet_perf] = train(tempNet,...
                inputs_train(:,cluster_i__ind{1,j}),...
                targ_train(:,cluster_i__ind{1,j}),...
                [],[],errorWeights);
            expertsNN{1,j} = tempNet;
            expertsNN{2,j} = tempNet_perf;
            Experts_perf_mat(j,i) = tempNet_perf.best_perf;
        else % give warning on the empty group
            % NOTE: is the group is empty, than the expert will not train!!
            warning(['at iter #',num2str(i),' the cluster of expert #',...
                num2str(j),' is empty']);
            expertsNN{2,j} = NaN;
            Experts_perf_mat(j,i) = NaN;
        end
    end
            
    % train the gate using f_h as targets:
    %       (minimize MSE between 'g' and 'f_h')
    [gateNet,gateNet_perf] = train(gateNet,inputs_train,fh);
    gateNN_perf_vec(1,i) = gateNet_perf.best_perf;
    
    % keep track on the test group Gate change 
    if i > 1
        [~,~,g_changes(i,:)] = ...
            obj.MoE_check_gate_change(gateOut_old,gateOut_test);
    end
    gateOut_old = gateOut_test;
    
    % plot the gate performance (MSE)
    if obj.disp_information
        plot(ax2,i,gateNet_perf.best_perf,'k','Marker','o');
    end
    
end

   
obj.MoE.expertsNN = expertsNN;
obj.MoE.gateNet = gateNet;
obj.MoE.gate_changes = g_changes; % rows-iteration num , col-mse chage per sample
obj.MoE.Moe_perf_over_iter = Moe_perf_over_iter;
obj.MoE.gateTraniData.gateNN_perf_vec = gateNN_perf_vec;
obj.MoE.out_from_train = netOut_train;
obj.MoE.out_from_valid = MoE_out_valid;

% test network to check test error:
    [MoE_out_test,~,~,~] = obj.MoE_testNet(inputs_test,expertsNN,...
    gateNet,type);
obj.MoE.out_from_test = MoE_out_test;
[Moe_MSE_on_test_temp,Rsquar] = ...
    obj.MoE_perf_calc(targ_test,MoE_out_test,0,0);

obj.MoE.expertsTrainData.expert_i_GroupSize = expert_i_GroupSize;
obj.MoE.expertsTrainData.Experts_perf_mat = Experts_perf_mat;
obj.MoE.expertsTrainData.emptyGroupIndecator = emptyGroupIndecator;

obj.MoE.Moe_RMSE_on_test = sqrt(Moe_MSE_on_test_temp);
obj.MoE.RsquarTest = Rsquar;

if obj.disp_information
    disp(['runTime of training: ',num2str(toc)]);
    disp([' MoE perf (MSE) on test group: ',num2str(Moe_MSE_on_test_temp)]);
end

end

