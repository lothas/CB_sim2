function [obj] = MoE_init(obj)
% this function initilize the MoE network.
% this function is being called by "MoE_train" function

%% initiate some stuff before
% set the experts colors for the graph colors 
switch obj.expertCount
    case {2,3} % in case of small number of expert, make colors clear:
        obj.colors = [1,0,0;0,1,0;0,0,1]; % RGB
    otherwise
        obj.colors = rand(obj.expertCount,3); % random colors
end

% defining the Legend with the experts names
legendNames_temp = cell(1,obj.expertCount);
for i=1:obj.expertCount
    legendNames_temp{1,i} = ['#',num2str(i),' expert'];
end
obj.legendNames = legendNames_temp;

%% continue to ini:
num_of_samples = size(obj.targets,2);

expertsNN = cell(1,obj.expertCount);

% define experts:
for j=1:obj.expertCount
    expertsNN{1,j} = feedforwardnet(obj.architecture);
    expertsNN{1,j}.trainParam.showWindow = 0; % dont show training window
    expertsNN{1,j}.trainParam.epochs = obj.MaxEpochs;
end

% initial training (to initilazied the NN)
% train each experts on a random sample of 10 points to initialze the
% weights.

if num_of_samples < 10
    error('num of train samples need to be larger than 10!');
end

for j=1:obj.expertCount
    initTrain = randsample(1:num_of_samples,10);
    [expertsNN{1,j}, ~] = train(expertsNN{1,j},...
        obj.inputs(:,initTrain),...
        obj.targets(:,initTrain));
end

% define and initialize gate Network:
g0 = softmax(rand(obj.expertCount,num_of_samples));
fh = obj.MoE_calc_fh(obj.inputs,obj.targets,expertsNN,g0);

initTrain = randsample(1:num_of_samples,10);
gateNet = trainSoftmaxLayer(obj.inputs(:,initTrain),fh(:,initTrain),...
    'LossFunction','mse','MaxEpochs',100,'ShowProgressWindow',false);

% gateNet = feedforwardnet(obj.my_MoE_out.GateHidLayer);
% layerNum = numel(obj.my_MoE_out.GateHidLayer); % number of layers
% % change the output process function to remap the output to be between 0to1
% gateNet.outputs{layerNum+1}.processParams{2}.ymin = 0;
% % change the output activation function to be 'softmax'
% gateNet.layers{layerNum+1}.transferFcn = 'softmax';
% gateNet.divideMode = 'none';
% gateNet.trainParam.showWindow = 0;
% % gateNet = train(gateNet,obj.sampl_train,fh);
% gateNet = train(gateNet,obj.sampl_train(:,1:10),mapminmax(fh(:,1:10),0,1)); % initilize faster on less samples?

obj.expertsNN = expertsNN;
obj.gateNet = gateNet;

end

