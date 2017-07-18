function [out] = NN_Perf_over_samples_num(sampl,targ,NumOfRepeats,HiddenN,...
    numOfSampl,train_or_plot,outputFileName)
% this function is calculating the NN performance over the number of
% neurons in the hidden layer.

% inputs:
% *) 'sampl' - the NN inputs
% *) 'targ' - the NN targets
% *) 'NumOfRepeats' - the number of times that we train each network (to
%                       get the statistical error).
% *) 'HiddenN' - the number of neurons in the hidden layer
% *) 'numOfSampl' - vector containing the number of samples to train
%                   ([10k,30k,...)
% *) 'train_or_plot' - #) 'train' - train the NN and save data to file
%					   #) 'plot' - plot the results on a graph
%					   #) 'text' - output the results in text


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

if max(numOfSampl) > samplNum*trainRatio
    error('not enougth samples...');
end

switch train_or_plot
    case 'train'
	
		disp('start the training...');
		
        netMseTrain = zeros(NumOfRepeats,length(HiddenN));
        netMseValidation = zeros(NumOfRepeats,length(HiddenN));
        netMseTest = zeros(NumOfRepeats,length(HiddenN));
        netTrainTime = zeros(NumOfRepeats,length(HiddenN));
        
        % Start Timer
        tic
        
        for i=1:length(numOfSampl)
		
			disp(['NN with ',num2str(numOfSampl(1,i)),' samples: ']);
			
            for j=1:NumOfRepeats

				train_start_time = toc;
				
				% Shuffle the samples (matrix rows):
				[trainInd,valInd,testInd] = dividerand(samplNum,trainRatio,valRatio,testRatio);
                
				
				% prepare the Neural Network:
                net = feedforwardnet(HiddenN);
                
				% Neural Network training parameters:
                net.trainParam.showWindow = false; % dont show training window
				net.divideFcn = 'divideind';
                net.divideParam.trainInd = trainInd(1:numOfSampl(1,i));
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
        out.samplesNum = numOfSampl;
        out.order_of_perf = {'1st col - train','2nd col - valid','3rd col - test'};
        
        % save the performances:
        out.netMseTrain = netMseTrain;
        out.netMseValidation = netMseValidation;
        out.netMseTest = netMseTest;
        
        meanMseTrain = mean(netMseTrain,1);
        meanMseValidation = mean(netMseValidation,1);
        meanMseTest = mean(netMseTest,1);
        out.NN_Mean_over_sampl_num = [meanMseTrain;meanMseValidation;meanMseTest];

        stdMseTrain = std(netMseTrain,0,1);
        stdMseValidation = std(netMseValidation,0,1);
        stdMseTest = std(netMseTest,0,1);
        out.NN_stdev_over_sampl_num = [stdMseTrain;stdMseValidation;stdMseTest];
        
        % save the training times:
        out.training_time_sec = netTrainTime;
        out.meanTrainTime = mean(netTrainTime,1);
        out.stdTrainTime = std(netTrainTime,0,1);

        save(outputFileName,'out')
        
    case 'plot'
        load(outputFileName,'out');
        
        graph_title = {'network perf_{(MSE)} over #samples',...
            ['NN with ',num2str(out.HiddenN),' hidden neurons']};

        
        meanMseTrain = out.NN_Mean_over_sampl_num(1,:);
        meanMseValidation = out.NN_Mean_over_sampl_num(2,:);
        meanMseTest = out.NN_Mean_over_sampl_num(3,:);
        
        stdMseTrain = out.NN_stdev_over_sampl_num(1,:);
        stdMseValidation = out.NN_stdev_over_sampl_num(2,:);
        stdMseTest = out.NN_stdev_over_sampl_num(3,:);
        
        figure;
        errorbar(out.samplesNum,meanMseTrain,stdMseTrain); hold on
        errorbar(out.samplesNum,meanMseValidation,stdMseValidation);
        errorbar(out.samplesNum,meanMseTest,stdMseTest); grid minor;
        legend('Train group','validation group','Test group');
        title(graph_title);
        xlabel('number of training samples');
        ylabel('MSE');
        set(gca,'FontSize',13);
		savefig('figure6_NN_perf_over_samplNum')
        
        figure;
        errorbar(out.samplesNum,out.meanTrainTime,out.stdTrainTime);
        title('train time over number of training samples'); grid minor;
        xlabel('number of training samples');
        ylabel('Time [sec]');
        set(gca,'FontSize',13);
        savefig('figure6b_trainTime_over_samplNum')
        
        out=[];
        
    case 'text'
        load(outputFileName,'out');
        
        numOfNNtypes = length(out.samplesNum);
        if numOfNNtypes == 1 
            % if we check only one type of NN, show all groups
            rowNames = {'train';'validation';'test'};
            meanMse = out.NN_Mean_over_sampl_num;
            stdMse = out.NN_stdev_over_sampl_num;
            resTable = table(meanMse,stdMse,'RowNames',rowNames)
        end
        if numOfNNtypes > 1
            % show only test perf
           disp('Results on test group:');
           rowNames = cell(numOfNNtypes,1);
           for i=1:numOfNNtypes
               rowNames{i,1} = [num2str(out.samplesNum(1,i)),' samples'];
           end
           meanMse = (out.NN_Mean_over_sampl_num(1,:))';
           stdMse = (out.NN_stdev_over_sampl_num(1,:))';
           resTable = table(meanMse,stdMse,'RowNames',rowNames)
        end
        
    otherwise
		warning('option are: 1)"train" 2)"plot" 3)"text"');
        error('invalid input "train_or_plot"');
end


end

