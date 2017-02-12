function [ NN_Mean_over_sampl_num,NN_stdev_over_sampl_num ] = NN_Perf_over_sampl_num( NumOfRepeats,HiddenN,dataPointsNum,NNinput,NNtarg,graphGO )
% this function is calculating the NN performance over the number of
% neurons in the hidden layer.

% inputs:
% 1) 'NumOfRepeats' - the amount of times that we train each network (to
%                       get the error bars).
% 2) 'HiddenN' - the number of neurons in the hidden layer
% 3) 'dataPointsNum' - the number of samples for the NN
% 4) 'NNinput' - the input matrix to NN
% 5) 'NNtarg' - the NN targets
% 6) 'graphGO' - to plot the error burs or not.

% ouputs:
% 1) 'NN_Mean_over_sampl_num' - the NN mean performance over the amount of
%                            data points.
% 2) 'NN_stdev_over_sampl_num' - the NN stdev of the performance over the amount of
%                            data points.

netMseTrain = zeros(NumOfRepeats,length(dataPointsNum));
netMseValidation = zeros(NumOfRepeats,length(dataPointsNum));
netMseTest = zeros(NumOfRepeats,length(dataPointsNum));

for i=1:length(dataPointsNum)
    for j=1:NumOfRepeats
        net = feedforwardnet(HiddenN);
        dataInd4Train = randsample(length(NNtarg),dataPointsNum(1,i));
        sampl4train = NNinput(:,dataInd4Train);
        targ4train = NNtarg(1,dataInd4Train);
        [net, tr] = train(net, sampl4train, targ4train);
        netMseTrain(j,i) = tr.best_perf;
        netMseValidation(j,i) = tr.best_vperf;
        netMseTest(j,i) = tr.best_tperf;
        
        clear net tr dataInd4Train sampl4train targ4train
    end
end

% save('NNdata.mat','netMseTrain','netMseValidation','netMseTest','HiddenN')
% save('NNdata_numPoints.mat','netMseTrain','netMseValidation','netMseTest','HiddenN','dataPointsNum')
meanMseTrain = mean(netMseTrain,1);
meanMseValidation = mean(netMseValidation,1);
meanMseTest = mean(netMseTest,1);
NN_Mean_over_sampl_num = [meanMseTrain,meanMseValidation,meanMseTest];

stdMseTrain = std(netMseTrain,0,1);
stdMseValidation = std(netMseValidation,0,1);
stdMseTest = std(netMseTest,0,1);
NN_stdev_over_sampl_num= [stdMseTrain,stdMseValidation,stdMseTest];

if graphGO
    figure;
    errorbar(dataPointsNum,meanMseTrain,stdMseTrain); hold on
    errorbar(dataPointsNum,meanMseValidation,stdMseValidation);
    errorbar(dataPointsNum,meanMseTest,stdMseTest);
    legend('Train group','validation group','Test group');
    title('network performance over samples num');
    xlabel('samples num');
    ylabel('MSE');
end

end


