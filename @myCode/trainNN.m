function [obj] = trainNN(obj,showWindow)
% this function trains a neural network using NNtoolbox

% TODO: add more control on the training parameters

% 'showWindow' - wether to show the training GUI or not 

obj.NN.net.trainParam.showWindow = showWindow; % dont show training window

disp('training neural network...');
[obj.NN.net, obj.NN.net_perf] = train(obj.NN.net, obj.sampl_train, obj.targ_train);

obj.NN.out_from_train = obj.NN.net(obj.sampl_train);
obj.NN.out_from_test = obj.NN.net(obj.sampl_test);

figure;
plotregression(obj.targ_train,obj.NN.out_from_train,'train',...
    obj.targ_test,obj.NN.out_from_test,'test');

% Calc test MSE:
[errMSE,~] = obj.NN_perf_calc(obj.targ_test,obj.NN.out_from_test,0,0,'test');
disp(['NN MSE on test group is: ',num2str(errMSE)]);

obj.NN.MSE_test_perf = errMSE;
end

