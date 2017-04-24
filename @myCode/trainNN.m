function [obj] = trainNN(obj,showWindow)
% this function trains a neural network using NNtoolbox

% TODO: add more control on the training parameters

% 'showWindow' - wether to show the training GUI or not 

obj.NN.net.trainParam.showWindow = showWindow; % dont show training window

trainSize = size(obj.sampl_train,2);
validSize = size(obj.sampl_valid,2);

sampl = horzcat(obj.sampl_train,obj.sampl_valid,obj.sampl_test);
targ = horzcat(obj.targ_train,obj.targ_valid,obj.targ_test);
obj.NN.net.divideFcn = 'divideind';
obj.NN.net.divideParam.trainInd = 1:trainSize;
obj.NN.net.divideParam.valInd   = trainSize+1:(trainSize+validSize);
obj.NN.net.divideParam.testInd  = (trainSize+validSize+1):size(sampl,2);
 
%  obj.NN.net.divideMode = 'none'; % all data to training
 
if obj.disp_information
    disp('training neural network...');
end

% Training:
%   also save a checkpoint file every 600 seconds (10 minutes). just in case
%   the computer will crash.
[obj.NN.net, obj.NN.net_perf] = train(obj.NN.net, sampl, targ,...
    'CheckpointFile','MyCheckpoint.mat','CheckpointDelay',600);

obj.NN.out_from_train = obj.NN.net(obj.sampl_train);
obj.NN.out_from_valid = obj.NN.net(obj.sampl_valid);
obj.NN.out_from_test = obj.NN.net(obj.sampl_test);

% Calc test MSE:
[errMSE,Rsquar] = obj.NN_perf_calc(obj.targ_test,obj.NN.out_from_test,0,0,'test');
obj.NN.MSE_test_perf = errMSE;
obj.NN.RsquarTest = Rsquar;

if obj.disp_information
    figure;
    plotregression(obj.targ_train,obj.NN.out_from_train,'train',...
        obj.targ_valid,obj.NN.out_from_valid,'validation',...
        obj.targ_test,obj.NN.out_from_test,'test');


    disp(['NN MSE on test group is: ',num2str(errMSE)]);

end
obj.NN.MSE_test_perf = errMSE;
end

