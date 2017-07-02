function obj = NN_multi_training(obj,varargin)
% this function trains a neural network using NNtoolbox

% 'showWindow' - wether to show the training GUI or not 

if obj.disp_information
    disp('training multi neural networks...');
    tic
end

switch nargin
    case 4
        obj.multi_NN.hiddenNeuronNum = varargin{1};
        hidN_vec = varargin{1};
        epochsNum = varargin{2};
        showWindow = varargin{3};
    case 3
        obj.multi_NN.hiddenNeuronNum = varargin{1};
        hidN_vec = varargin{1};
        epochsNum = varargin{2};
        showWindow = false;
    otherwise
        disp('invalid number of inputs...')
        error('the inputs order is: 1)hidNeurons, 2)numEphocs, 3)showWind');
end

%% get inputs and targets
sampl = obj.sampl;
targ = obj.targ;

%% prepare for training
outputsNum = size(targ,1);
NNs = cell(1,outputsNum);
NNs_trainPerf = cell(1,outputsNum);
NNs_testPerf = zeros(1,outputsNum);

for j=1:outputsNum
    tempNet = feedforwardnet(hidN_vec);
    tempNet.trainParam.epochs = epochsNum;
    tempNet.trainParam.showWindow = showWindow; % dont show training window
    
    [tempNet,temp_tr] = train(tempNet, sampl, targ(j,:));
    
    NNs_testPerf(1,j) = temp_tr.best_tperf;
    NNs{1,j} = tempNet;
    NNs_trainPerf{1,j} = temp_tr;
end

obj.multi_NN.NNs = NNs;
obj.multi_NN.net_perf = NNs_trainPerf;
obj.multi_NN.outNames = obj.outputNames;

if obj.disp_information
    % get the number of weights:
    obj.multi_NN.num_of_weights = tempNet.numWeightElements;
    disp(['training time was: ',num2str(toc), '[sec]']);
    disp(['the number of samples is: ',num2str(size(sampl,2))]);
    disp(['the number of weights in each NN: ',...
        num2str(obj.multi_NN.num_of_weights)]);
    disp(['the ratio between #sampl/#weights =  ',...
        num2str(size(sampl,2)/obj.multi_NN.num_of_weights)]);
    
    for j=1:outputsNum
       disp(['MSE test perf for output: ',obj.outputNames{1,j},...
           ' is ',num2str(NNs_testPerf(1,j))]);
    end
    
%     NN_multi_perf_plots(obj)
end

end

