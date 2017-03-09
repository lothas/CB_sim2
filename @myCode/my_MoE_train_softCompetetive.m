function [obj] = my_MoE_train_softCompetetive(obj)
%this function train a "Mixture of Experts" newural networks
% % NOTE: same for the "my_MoE_train.m" but only for competetiveFlag=2 and
% with better visualization.

numToShow = 100; % the #samples to show in the online regression plot.

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

sampl_train = obj.sampl_train;
targ_train = obj.targ_train;
sampl_valid = obj.sampl_valid;
targ_valid = obj.targ_valid;
sampl_test = obj.sampl_test;
targ_test = obj.targ_test;

% take only few points for the regression plot:
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
errMat = zeros(expertCount,num_of_train_samples); % error matrix (each row - targets)
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
    plot(ax1,i,MSE_train,'b','Marker','o');
    plot(ax1,i,MSE_valid,'g','Marker','o');
    plot(ax1,i,MSE_test,'r','Marker','o');
    
    % plot and update the regression plot:
    [g_max,g_max_ind] = max(gateOut_test,[],1);
    ouputs = MoE_out_test;
    targets = targ_test;
    cla(ax3);
    for k=1:obj.expertCount
        for n=1:numToShow
            line([0,1],[0,1],'color',[0,0,0],'linewidth',2,'Parent',ax3);
            ind = randSampl_ind(1,n);
            if g_max_ind(1,ind) == k
                if g_max(1,ind) > 0.5
                    plot(ax3,targets(1,ind),ouputs(1,ind),'k-o','MarkerFaceColor',obj.colors(k,:));
                else
                    plot(ax3,targets(1,ind),ouputs(1,ind),'k-o');
                end
            end
        end
    end
    
    cla(ax4);
    bplot = bar(ax4,(gateOut_test(:,randSampl_ind))','stacked');
    for k=1:length(randSampl_ind)
      set(bplot(k),'facecolor',obj.colors(k,:))
    end
    AX=legend(H, {'a','b','c','d','e','f'}, 'Location','Best','FontSize',8);
    
    % calc MoE performance to check wether to stop the training:
    [Moe_perf_over_iter(1,i),~] = obj.NN_perf_calc(targ_valid,MoE_out_valid,0,0);
    if Moe_perf_over_iter(1,i) < 0.0000001 % stopping condition on error
        disp('reached below the desired error');
        break;
    end
    if (i > 11) && (mean(Moe_perf_over_iter(1,(i-10):i)) < 0.00000001) % stopping condition on error gradient
        disp('reached below the desired error gradient');
        break;
        % TODO: use 'diff' to calculate the gradient.
    end
    
   % check the cluster size of each expert: 
   for j=1:expertCount % check the size of each cluster
        expert_i_GroupSize(j,i) = length(cluster_i__ind{1,j}); 
   end
   
   % run each expert on the entire data and calc "fh":
   g = gateOut;
   fh = obj.calc_fh(expertsNN,g);
   
    % training the experts:
    %   train each expert only on the samples which belongs to it.
    %      also train the experts with weights given to samples by f_h 
    for j=1:expertCount
        tempNet = expertsNN{1,j};
        errorWeights_all0 = {fh(j,:)};
        errorWeights_all = errorWeights_all0{1,1};
        errorWeights = errorWeights_all(1,cluster_i__ind{1,j});
        
        [tempNet, tempNet_perf] = train(tempNet,...
            sampl_train(:,cluster_i__ind{1,j}),...
            targ_train(:,cluster_i__ind{1,j}),...
            [],[],errorWeights);
            % with soft competetive the expert's cluster can't be empty
            % (it's not likely)
        expertsNN{1,j} = tempNet;
        expertsNN{2,j} = tempNet_perf;
        Experts_perf_mat(j,i) = tempNet_perf.best_perf;
    end
            

    % train the gate using f_h as targets:
    %       (minimize MSE between 'g' and 'f_h')
    [gateNet,gateNet_perf] = train(gateNet,sampl_train,fh);
    plot(ax2,i,gateNet_perf.best_perf,'k','Marker','o');
    gateNN_perf_vec(1,i) = gateNet_perf.best_perf;
    
%     % % % ploting the expert's weights:
%     ax_temp = gca;
%     title_temp = (['expert #',num2str(1)]);
%     axes(ax2);
%     obj.NN_weights_matrix_visualize(expertsNN{1,1},title_temp);
%     axes(ax_temp);
    
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

