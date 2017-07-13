function [out] = NN_Perf_over_Hid_Neuron_Num(sampl,targ,NumOfRepeats,HiddenN,train_or_plot,...
    keepRatioConstant,outputFileName)
% this function is calculating the NN performance over the number of
% neurons in the hidden layer.

% inputs:
% *) 'sampl' - the NN inputs
% *) 'targ' - the NN targets
% *) 'NumOfRepeats' - the number of times that we train each network (to
%                       get the statistical error).
% *) 'HiddenN' - the number of neurons in the hidden layer
% *) 'train_or_plot' - #) 'train' - train the NN and save data to file
%					   #) 'plot' - plot the results on a graph
%					   #) 'text' - output the results in text
% *) 'keepRatioConstant' - if 'true'- keep the ratio between the number of
%                               neurons to the number of samples constant.
%                          if 'false' - take the number of samples from the
%                                       class (no change)

% ouputs:
% 1) 'NN_Mean_over_HN_num' - the NN mean performance over the amount of
%                            neurons in the hidden layer.
% 2) 'NN_stdev_over_HN_num' - the NN stdev of the performance over the amount of
%                            neurons in the hidden layer.


samplNum = size(sampl,2); % the number of samples in total
num_of_inputs = size(sampl,1);
num_of_outputs = size(targ,1);

trainRatio = 0.7; 
valRatio = 0.15;
testRatio = 1 - trainRatio - valRatio;

switch train_or_plot
    case 'train'
	
		disp('start the training...');
		
        netMseTrain = zeros(NumOfRepeats,length(HiddenN));
        netMseValidation = zeros(NumOfRepeats,length(HiddenN));
        netMseTest = zeros(NumOfRepeats,length(HiddenN));
        netTrainTime = zeros(NumOfRepeats,length(HiddenN));
        
        % Start Timer
        tic
        
        for i=1:length(HiddenN)
		
			disp(['NN with ',num2str(HiddenN(1,i)),' hidden neurons: ']);
			
            for j=1:NumOfRepeats

				train_start_time = toc;
				
				% Shuffle the samples (matrix rows):
				[trainInd,valInd,testInd] = dividerand(samplNum,trainRatio,valRatio,testRatio);
				train_size = length(trainInd);
				
                
                % calc the number of weights:
                num_of_weights = ( (num_of_inputs+1) * HiddenN(1,i) )  +...
                    ( (HiddenN(1,i)+1) * num_of_outputs ); 

                % If to take the training samples as is, or to reduce their
                %       number to only 500 times more than the NN weights:
                ratio = 500; %times more samples than weights
                % calc if we have enough samples to reduce the number
                enough_sampl = train_size > (ratio*num_of_weights);
                if keepRatioConstant && enough_sampl
					% take only the 'n' first sampl from the hidden group
                    trainInd = trainInd(1:floor(ratio*num_of_weights));
                end
                
				if ~enough_sampl
					% alert if we don't a the desired amount of training samples
					warning('not enough samples to reduce... took everything we have');
				end
				
				% prepare the Neural Network:
                net = feedforwardnet(HiddenN(1,i));
                
				% Neural Network training parameters:
                net.trainParam.showWindow = false; % dont show training window
				net.divideFcn = 'divideind';
                net.divideParam.trainInd = trainInd;
                net.divideParam.valInd   = valInd;
                net.divideParam.testInd  = testInd;
                
                % NN training
                [~, tr] = train(net, sampl, targ);
                netMseTrain(j,i) = tr.best_perf;
                netMseValidation(j,i) = tr.best_vperf;
                netMseTest(j,i) = tr.best_tperf;
                
                clear tr
				
				train_end_time = toc;
				train_time = train_end_time - train_start_time;
                
                netTrainTime(j,i) = train_time;
				
				disp(['    ',num2str(j),' out of ',num2str(NumOfRepeats),...
					' train time = ',num2str(train_time),'[sec]']);
				
            end
        end

        out.HiddenN = HiddenN;
        out.samplesNum = samplNum;
        out.order_of_perf = {'1st col - train','2nd col - valid','3rd col - test'};
        
        % save the performances:
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
        
        % save the training times:
        out.training_time_sec = netTrainTime;
        out.meanTrainTime = mean(netTrainTime,1);
        out.stdTrainTime = std(netTrainTime,0,1);
        
        save(outputFileName,'out')
        
    case 'plot'
        load(outputFileName,'out');
        
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
        errorbar(out.HiddenN,meanMseTest,stdMseTest); grid minor;
        legend('Train group','validation group','Test group');
        title(graph_title);
        xlabel('Hidden Neuron Num');
        ylabel('MSE');
        set(gca,'FontSize',13);
		
        figure;
        errorbar(out.HiddenN,out.meanTrainTime,out.stdTrainTime); grid minor;
        title('train time over number of hidden neurons');
        xlabel('number of neurons');
        ylabel('Time [sec]');
        set(gca,'FontSize',13);
        
        out=[];
        
    case 'text'
        load(outputFileName,'out');
        
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
		warning('option are: 1)"train" 2)"plot" 3)"text"');
        error('invalid input "train_or_plot"');
end


end

