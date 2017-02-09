function [out] = NN_Perf_over_sampl_num( NumOfRepeats,HiddenN,dataPointsNum,NNinput,NNtarg,train_or_plot )
% this function is calculating the NN performance over the number of
% neurons in the hidden layer.

% inputs:
% 1) 'NumOfRepeats' - the amount of times that we train each network (to
%                       get the error bars).
% 2) 'HiddenN' - the number of neurons in the hidden layer
% 3) 'dataPointsNum' - the number of samples for the NN
% 4) 'NNinput' - the input matrix to NN
% 5) 'NNtarg' - the NN targets

switch train_or_plot
    case 'train'
        netMseTrain = zeros(NumOfRepeats,length(dataPointsNum));
        netMseValidation = zeros(NumOfRepeats,length(dataPointsNum));
        netMseTest = zeros(NumOfRepeats,length(dataPointsNum));

        for i=1:length(dataPointsNum)
            for j=1:NumOfRepeats
                net = feedforwardnet(HiddenN);
                dataInd4Train = randsample(size(NNtarg,2),dataPointsNum(1,i));
                sampl4train = NNinput(:,dataInd4Train);
                targ4train = NNtarg(:,dataInd4Train);
                [~, tr] = train(net, sampl4train, targ4train);
                netMseTrain(j,i) = tr.best_perf;
                netMseValidation(j,i) = tr.best_vperf;
                netMseTest(j,i) = tr.best_tperf;

                clear tr dataInd4Train sampl4train targ4train
            end
        end
        
        out.HiddenN = HiddenN;
        out.samplesNum = dataPointsNum;
        
        out.netMseTrain = netMseTrain;
        out.netMseValidation = netMseValidation;
        out.netMseTest = netMseTest;

        meanMseTrain = mean(netMseTrain,1);
        meanMseValidation = mean(netMseValidation,1);
        meanMseTest = mean(netMseTest,1);
        out.NN_Mean_over_sampl_num = [meanMseTrain,meanMseValidation,meanMseTest];

        stdMseTrain = std(netMseTrain,0,1);
        stdMseValidation = std(netMseValidation,0,1);
        stdMseTest = std(netMseTest,0,1);
        out.NN_stdev_over_sampl_num= [stdMseTrain,stdMseValidation,stdMseTest];
        
        out.order_of_perf = {'1st col - train','2nd col - vald','3rd col - test'};
        
        save('NN_Perf_over_sampl_num_result.mat','out')
        
    case 'plot'
        load('NN_Perf_over_sampl_num_result.mat','out');
        
        figure;
        errorbar(out.samplesNum,out.meanMseTrain,out.stdMseTrain); hold on
        errorbar(out.samplesNum,out.meanMseValidation,out.stdMseValidation);
        errorbar(out.samplesNum,out.meanMseTest,out.stdMseTest);
        legend('Train group','validation group','Test group');
        title('network performance over samples num');
        xlabel('samples num');
        ylabel('MSE');
        
    otherwise
        error('invalid input "train_or_plot"');
end


end


