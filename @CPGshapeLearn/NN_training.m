function obj = NN_training(obj,varargin)
% this function trains a neural network using NNtoolbox

% TODO: add more control on the training parameters

% 'showWindow' - wether to show the training GUI or not 

switch nargin
    case 5
        obj.NN.hiddenNeuronNum = varargin{1};
        obj.NN.net = feedforwardnet(varargin{1});
        obj.NN.net.trainParam.epochs = varargin{2};
        showWindow = varargin{3};
        useGPU_flag = varargin{4};
    case 3
        obj.NN.hiddenNeuronNum = varargin{1};
        obj.NN.net = feedforwardnet(varargin{1});
        obj.NN.net.trainParam.epochs = varargin{2};
        showWindow = false;
        useGPU_flag = false;
    otherwise
        disp('invalid number of inputs...')
        error('the inputs order is: 1)hidNeurons, 2)numEphocs, 3)showWind, 4)useGPU');
end

obj.NN.usingGPU = useGPU_flag;
obj.NN.net.trainParam.showWindow = showWindow; % dont show training window

sampl = obj.sampl;
targ = obj.targ;

if obj.disp_information
    disp('training neural network...');
    tic
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

% get the number of weights:
obj.NN.num_of_weights = obj.NN.net.numWeightElements;

% Calc test MSE for each targ(output):
testInd = obj.NN.net_perf.testInd;
targ_test = targ(:,testInd);
output_test = obj.NN.net(sampl(:,testInd));
numOfOut = size(targ,1);
MSE_test = zeros(numOfOut,1);
for j=1:numOfOut
    [MSE_test(j,1),~] = obj.NN_perf(targ_test(j,:),output_test(j,:));
end
obj.NN.MSE_test_perf = MSE_test;

if obj.disp_information
    disp(['training time was: ',num2str(toc), '[sec]']);
    disp(['the number of samples is: ',num2str(size(sampl,2))]);
    disp(['the number of network weights: ',num2str(obj.NN.net.numWeightElements)]);
    disp(['the ratio between #sampl/#weights =  ',...
        num2str(size(sampl,2)/obj.NN.net.numWeightElements)]);
   
    disp(' ');
    errs = obj.NN.MSE_test_perf;
    disp('the MSE err for each output:');
    for i=1:length(obj.outputNames)
        disp([obj.outputNames{1,i},' = ',num2str(errs(i))]);
    end
end

end

