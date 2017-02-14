function [obj] = My_MoE_init(obj)
% this function initilize the MoE network.
% this function is being called by "My_MoE_train" function

if isempty(obj.sampl_train) || isempty(obj.targ_train)
    error('data was not splited to train and test groups');
end

num_of_train_samples = size(obj.sampl_train,2);

expertsNN = cell(2,obj.expertCount);

% define experts:
for j=1:obj.expertCount
    expertsNN{1,j} = feedforwardnet(obj.my_MoE_out.ExpertHidLayer);
    expertsNN{1,j}.trainParam.showWindow = 0; % dont show training window
    expertsNN{1,j}.trainParam.epochs = obj.my_MoE_out.maxEphocs;
    expertsNN{1,j}.divideMode = 'none'; % all data to training
end

% initial training (to initilazied the NN)
% train each experts on a random sample of 500 points to initialze the
% weights.

if num_of_train_samples < 500
    error('num of train samples need to be larger than 500!');
end

for j=1:obj.expertCount
    initTrain = randsample(1:num_of_train_samples,500);
    [expertsNN{1,j}, expertsNN{2,j}] = train(expertsNN{1,j}, obj.sampl_train(:,initTrain), obj.targ_train(:,initTrain));
end

% define and initialize gate Network:
for j=1:obj.expertCount % run the data throught the experts to get initial clustering
    tempNet = expertsNN{1,j};
    outMat(j,:) = tempNet(obj.sampl_train);
    errMat(j,:) = outMat(j,:) - obj.targ_train;
end
seMat = errMat.^2; % squar error

% [~,best_expert_ind] = min(seMat,[],1);
% for j=1:obj.expertCount % find targets based on experts perf
%     gateNet_targ(j,:) = (best_expert_ind == j);
% end
% gateNet = patternnet(obj.my_MoE_out.GateHidLayer);
% gateNet.performParam.regularization = 0;
% gateNet.trainParam.showWindow = 0;
% gateNet = train(gateNet,obj.sampl_train,gateNet_targ);

g0 = rand(obj.expertCount,num_of_train_samples);
fh = zeros(obj.expertCount,num_of_train_samples);
for k=1:num_of_train_samples
    fh(:,k) = g0(:,k) .* exp(-0.5 .* seMat(:,k) );
    fh(:,k) = fh(:,k) ./ sum(fh(:,k),1);
end

gateNet = feedforwardnet(obj.my_MoE_out.GateHidLayer);
layerNum = numel(obj.my_MoE_out.GateHidLayer); % number of layers
% change the output process function to remap the output to be between 0to1
gateNet.outputs{layerNum+1}.processParams{2}.ymin = 0;
% change the output activation function to be 'softmax'
gateNet.layers{gateNet.numLayers}.transferFcn = 'softmax';
gateNet.divideMode = 'none';
gateNet.trainParam.showWindow = 0;
gateNet = train(gateNet,obj.sampl_train,fh);

obj.my_MoE_out.expertsNN = expertsNN;
obj.my_MoE_out.gateNet = gateNet;

end

