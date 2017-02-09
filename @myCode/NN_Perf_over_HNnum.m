function [out] = NN_Perf_over_HNnum( NumOfRepeats,HiddenN,NNinput,NNtarg,train_or_plot )
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

switch train_or_plot
    case 'train'
        netMseTrain = zeros(NumOfRepeats,length(HiddenN));
        netMseValidation = zeros(NumOfRepeats,length(HiddenN));
        netMseTest = zeros(NumOfRepeats,length(HiddenN));

        for i=1:length(HiddenN)
            for j=1:NumOfRepeats
                net = feedforwardnet(HiddenN(1,i));
                [~, tr] = train(net, NNinput, NNtarg);
                netMseTrain(j,i) = tr.best_perf;
                netMseValidation(j,i) = tr.best_vperf;
                netMseTest(j,i) = tr.best_tperf;

                clear tr
            end
        end

        out.HiddenN = HiddenN;
        out.samplesNum = size(NNtarg,2);

        out.netMseTrain = netMseTrain;
        out.netMseValidation = netMseValidation;
        out.netMseTest = netMseTest;
        
        meanMseTrain = mean(netMseTrain,1);
        meanMseValidation = mean(netMseValidation,1);
        meanMseTest = mean(netMseTest,1);
        out.NN_Mean_over_HN_num = [meanMseTrain,meanMseValidation,meanMseTest];

        stdMseTrain = std(netMseTrain,0,1);
        stdMseValidation = std(netMseValidation,0,1);
        stdMseTest = std(netMseTest,0,1);
        out.NN_stdev_over_HN_num = [stdMseTrain,stdMseValidation,stdMseTest];
        
        out.order_of_perf = {'1st col - train','2nd col - vald','3rd col - test'};
        
        save('NN_Perf_over_HNnum.mat','out')
        
    case 'plot'
        load('NN_Perf_over_HNnum.mat','out');
        
        figure;
        errorbar(out.HiddenN,out.meanMseTrain,out.stdMseTrain); hold on
        errorbar(out.HiddenN,out.meanMseValidation,out.stdMseValidation);
        errorbar(out.HiddenN,out.meanMseTest,out.stdMseTest);
        legend('Train group','validation group','Test group');
        title('network performance (MSE) over hidden neurons num (target=frequency)');
        xlabel('Hidden Neuron Num');
        ylabel('MSE');
        
        out=[];
    otherwise
        error('invalid input "train_or_plot"');
end


end

