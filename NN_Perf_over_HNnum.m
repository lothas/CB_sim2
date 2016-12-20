function [ NN_Mean_over_HN_num,NN_stdev_over_HN_num ] = NN_Perf_over_HNnum( NumOfRepeats,HiddenN,NNinput,NNtarg,graphGO )
% this function is calculating the NN performance over the number of
% neurons in the hidden layer.

% inputs:
% 1) 'NumOfRepeats' - the amount of times that we train each network (to
%                       get the error bars).
% 2) 'HiddenN' - the number of neurons in the hidden layer
% 3) 'NNinput' - the input matrix to NN
% 4) 'NNtarg' - the NN targets
% 5) 'graphGO' - to plot the error burs or not.

% ouputs:
% 1) 'NN_Mean_over_HN_num' - the NN mean performance over the amount of
%                            neurons in the hidden layer.
% 2) 'NN_stdev_over_HN_num' - the NN stdev of the performance over the amount of
%                            neurons in the hidden layer.

netMseTrain = zeros(NumOfRepeats,length(HiddenN));
netMseValidation = zeros(NumOfRepeats,length(HiddenN));
netMseTest = zeros(NumOfRepeats,length(HiddenN));

for i=1:length(HiddenN)
    for j=1:NumOfRepeats
        net = feedforwardnet(HiddenN(1,i));
        [net, tr] = train(net, NNinput, NNtarg);
        netMseTrain(j,i) = tr.best_perf;
        netMseValidation(j,i) = tr.best_vperf;
        netMseTest(j,i) = tr.best_tperf;
        
        clear net tr
    end
end

% save('NNdata.mat','netMseTrain','netMseValidation','netMseTest','HiddenN')
meanMseTrain = mean(netMseTrain,1);
meanMseValidation = mean(netMseValidation,1);
meanMseTest = mean(netMseTest,1);
NN_Mean_over_HN_num = [meanMseTrain,meanMseValidation,meanMseTest];

stdMseTrain = std(netMseTrain,0,1);
stdMseValidation = std(netMseValidation,0,1);
stdMseTest = std(netMseTest,0,1);
NN_stdev_over_HN_num = [stdMseTrain,stdMseValidation,stdMseTest];

if graphGO
    figure;
    errorbar(HiddenN,meanMseTrain,stdMseTrain); hold on
    errorbar(HiddenN,meanMseValidation,stdMseValidation);
    errorbar(HiddenN,meanMseTest,stdMseTest);
    legend('Train group','validation group','Test group');
    title('network performance (MSE) over hidden neurons num (target=frequency)');
    xlabel('Hidden Neuron Num');
    ylabel('MSE');
end

end

