function [obj] = trainNN(obj,varargin)
% this function trains a neural network using NNtoolbox

% TODO: add more control on the training parameters

% 'showWindow' - wether to show the training GUI or not 

switch nargin
    case 3
        showWindow = varargin{1};
        useGPU_flag = varargin{2};
    case 2
        showWindow = varargin{1};
        useGPU_flag = false;
    case 1
        showWindow = false;
        useGPU_flag = false;
    otherwise
        error('invalid number of inputs')
end

obj.NN.usingGPU = useGPU_flag;

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
 
% calc the number of weights:
obj.NN.num_of_weights = length(getwb(obj.NN.net));

if obj.disp_information
    disp('training neural network...');
end

if ~useGPU_flag
    % % % Training (NO GPU): % % %
    %   also save a checkpoint file every 1800 seconds (30 minutes). just in case
    %   the computer will crash.
%     [obj.NN.net, obj.NN.net_perf] = train(obj.NN.net, sampl, targ,...
%         'CheckpointFile','MyCheckpoint.mat','CheckpointDelay',1800);
    [obj.NN.net, obj.NN.net_perf] = train(obj.NN.net, sampl, targ);
else
    % % % Training (WITH GPU): % % %
    % note: If NET has the default training function trainlm,
    %   you see a warning that GPU calculations do not support Jacobian training,
    %   only gradient training. So the training function is automatically
    %   changed to the gradient training function trainscg.
    %   To avoid the notice, you can specify the function before training:
    %      "net.trainFcn = 'trainscg';" 
    [obj.NN.net, obj.NN.net_perf] = train(obj.NN.net, sampl, targ,...
        'useGPU','yes');
end
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

