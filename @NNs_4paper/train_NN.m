function obj = train_NN(obj,architecture,inputs,targets)
% create and train the NN

% 'architecture' - row vector with the NN architecture
% 'inputs' - NN inputs
% 'targets' - NN targets

targets_num = size(targets,1);

% Create and train the NN
net = feedforwardnet(architecture);
[net, ~] = train(net, inputs, targets);

% saving the 'net' object:
obj.NN.net = net;

% calculating NN perf for each target:
NNoutput = net(inputs);
perf = cell(1,targets_num);
for i=1:targets_num
    perf{1,i} = immse(NNoutput(i,:),targets(i,:));
end
obj.NN.NN_errMSE = perf;


% % % R^2 = 1 - (error variance/input variance)
% calc R^2
% err = Targets-NNoutput;
% errVar = var(err,0,2);
% inputVar = var(Targets,0,2);
% R_squar = 1-(errVar/inputVar);

end

