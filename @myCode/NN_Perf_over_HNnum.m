function [out] = NN_Perf_over_HNnum(obj,NumOfRepeats,HiddenN,train_or_plot,...
    keepRatioConstant)
% this function is calculating the NN performance over the number of
% neurons in the hidden layer.

% inputs:
% 1) 'NumOfRepeats' - the amount of times that we train each network (to
%                       get the error bars).
% 2) 'HiddenN' - the number of neurons in the hidden layer
% 3) 'train_or_plot' - to plot the error burs or not.
% 4) 'keepRatioConstant' - if 'true'- keep the ratio between the number of
%                               neurons to the number of samples constant.
%                          if 'false' - take the number of samples from the
%                                       class (no change)

% ouputs:
% 1) 'NN_Mean_over_HN_num' - the NN mean performance over the amount of
%                            neurons in the hidden layer.
% 2) 'NN_stdev_over_HN_num' - the NN stdev of the performance over the amount of
%                            neurons in the hidden layer.

switch train_or_plot
    case 'train'
        netMseTrain = zeros(NumOfRepeats,length(HiddenN));
        netMseValidation = zeros(NumOfRepeats,length(HiddenN));
        netMseTest = zeros(NumOfRepeats,length(HiddenN));
                
        for i=1:length(HiddenN)
            for j=1:NumOfRepeats
                % NOTE: this code shuffles the train and validation group in each
                % iteration. but keeps the test group fixed.

                sampl_train = obj.sampl_train;
                targ_train = obj.targ_train;
                
                % calc the number of weights:
                num_of_inputs = size(sampl_train,1);
                num_of_outputs = size(targ_train,1);
                num_of_weights = ( (num_of_inputs+1) * HiddenN(1,i) )  +...
                    ( (HiddenN(1,i)+1) * num_of_outputs ); 
                
                % If to take the training samples as is, or to reduce their
                %       number to only 100 times more than the NN weights:
                ratio = 100; %100 times more samples than weights
                % calc if we have enough samples to reduce the number
                enough_cond = size(sampl_train,2) > (ratio*num_of_weights);
                if keepRatioConstant && enough_cond
                    rand_ind = randsample(1:size(obj.sampl_train,2),...
                        ratio*num_of_weights);
                    sampl_train = sampl_train(:,rand_ind);
                    targ_train = targ_train(:,rand_ind);
                end
                
                trainSize = size(sampl_train,2);
                validSize = size(obj.sampl_valid,2);
                
                net = feedforwardnet(HiddenN(1,i));
                
                net.trainParam.showWindow = false; % dont show training window
                sampl = horzcat(sampl_train,obj.sampl_valid,obj.sampl_test);
                targ = horzcat(targ_train,obj.targ_valid,obj.targ_test);

                % set Train, Valid and Test groups:
                net.divideFcn = 'divideind';
                net.divideParam.trainInd = 1:trainSize;
                net.divideParam.valInd   = trainSize+1:(trainSize+validSize);
                net.divideParam.testInd  = (trainSize+validSize+1):size(sampl,2);
                
                % NN training
                [~, tr] = train(net, sampl, targ);
                netMseTrain(j,i) = tr.best_perf;
                netMseValidation(j,i) = tr.best_vperf;
                netMseTest(j,i) = tr.best_tperf;

                % %  shuffle for next time:
                % only for 2neuron symmetric case-
                obj = obj.shuffle_samples('onlyTrainAndValid'); 
                % for the other cases- dont keep test group the same-
%                 obj = obj.shuffle_samples('completeShuffle');
                
                clear tr
            end
        end

        out.HiddenN = HiddenN;
        out.samplesNum = size(targ,2);

        out.netMseTrain = netMseTrain;
        out.netMseValidation = netMseValidation;
        out.netMseTest = netMseTest;
        
        meanMseTrain = mean(netMseTrain,1);
        meanMseValidation = mean(netMseValidation,1);
        meanMseTest = mean(netMseTest,1);
        out.NN_Mean_over_HN_num = [meanMseTrain;meanMseValidation;meanMseTest];

        stdMseTrain = std(netMseTrain,0,1);
        stdMseValidation = std(netMseValidation,0,1);
        stdMseTest = std(netMseTest,0,1);
        out.NN_stdev_over_HN_num = [stdMseTrain;stdMseValidation;stdMseTest];
        
        out.order_of_perf = {'1st col - train','2nd col - vald','3rd col - test'};
        
        save('NN_Perf_over_HNnum.mat','out')
        
    case 'plot'
        load('NN_Perf_over_HNnum.mat','out');
        
        if keepRatioConstant
            % TODO: make sure that we have here enough samples to reduce.
            graph_title = {'network perf_{(MSE)} over #hidden_{N}',...
                '#train_{samples} = K * #Weights'};
        else
            graph_title = {'network perf(MSE) over #hidden_{neurons}',...
                ['#train_{samples} = ',num2str(size(obj.sampl_train,2))]};
        end
        
        meanMseTrain = out.NN_Mean_over_HN_num(1,:);
        meanMseValidation = out.NN_Mean_over_HN_num(2,:);
        meanMseTest = out.NN_Mean_over_HN_num(3,:);
        
        stdMseTrain = out.NN_stdev_over_HN_num(1,:);
        stdMseValidation = out.NN_stdev_over_HN_num(2,:);
        stdMseTest = out.NN_stdev_over_HN_num(3,:);
        
        figure;
        errorbar(out.HiddenN,meanMseTrain,stdMseTrain); hold on
        errorbar(out.HiddenN,meanMseValidation,stdMseValidation);
        errorbar(out.HiddenN,meanMseTest,stdMseTest);
        legend('Train group','validation group','Test group');
        title(graph_title);
        xlabel('Hidden Neuron Num');
        ylabel('MSE');
        
        out=[];
        
    case 'text'
        load('NN_Perf_over_HNnum.mat','out');
        
        numOfNNtypes = length(out.HiddenN);
        if numOfNNtypes == 1 
            % if we check only one type of NN, show all groups
            rowNames = {'train';'validation';'test'};
            meanMse = out.NN_Mean_over_HN_num;
            stdMse = out.NN_stdev_over_HN_num;
            resTable = table(meanMse,stdMse,'RowNames',rowNames)
        end
        if numOfNNtypes > 1
            % show only test perf
           disp('Results on test group:');
           rowNames = cell(numOfNNtypes,1);
           for i=1:numOfNNtypes
               rowNames{i,1} = [num2str(out.HiddenN(1,i)),' hidden neurons'];
           end
           meanMse = (out.NN_Mean_over_HN_num(1,:))';
           stdMse = (out.NN_stdev_over_HN_num(1,:))';
           resTable = table(meanMse,stdMse,'RowNames',rowNames)
        end
        
    otherwise
        error('invalid input "train_or_plot"');
end


end

