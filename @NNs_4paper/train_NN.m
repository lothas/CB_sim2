function obj = train_NN(obj,architecture,inputs,targets)
% create and train the NN

% 'architecture' - row vector with the NN architecture
% 'inputs' - NN inputs
% 'targets' - NN targets

targets_num = size(targets,1);

% Create and train the NN
net = feedforwardnet(architecture);

% net.trainFcn = 'trainlm';
% net.trainFcn = 'trainbr';
% net.trainFcn = 'trainscg';
% net.trainFcn = 'trainrp';

if false % train only on some of the samples (not all of them)
    ind = randsample(size(targets,2),50000);
    inputs = inputs(:,ind);
    targets = targets(:,ind);
end

[net, ~] = train(net, inputs, targets);

% saving the 'net' object:
obj.NN.net = net;

% calculating NN perf for each target:
NNoutput = net(inputs);
perf = cell(1,targets_num);
for i=1:targets_num
    perf{1,i} = sqrt(immse(NNoutput(i,:),targets(i,:)));
end
obj.NN.NN_errRMSE = perf;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xedges = linspace(0,0.25,100);
% Yedges = linspace(0,10,100);
% 
% figure; hold on;
% % scatter(tau_out_on_train_set,b_out_on_train_set,'.k');
% histogram2(NNoutput(1,:),NNoutput(2,:),Xedges,Yedges,...
%     'DisplayStyle','tile','ShowEmptyBins','on',...
%     'Normalization','pdf');
% % scatter(tau_out_on_train_set,b_out_on_train_set,'.k');
% title(sprintf('2D distribution of tau and b  \n on training data \n NN output'));
% xlabel('tau');
% ylabel('b');
% 
% figure; hold on;
% % scatter(tau_out_on_train_set,b_out_on_train_set,'.k');
% histogram2(targets(1,:),targets(2,:),Xedges,Yedges,...
%     'DisplayStyle','tile','ShowEmptyBins','on',...
%     'Normalization','pdf');
% % scatter(tau_out_on_train_set,b_out_on_train_set,'.k');
% title(sprintf('2D distribution of tau and b  \n on training data \n NN targets'));
% xlabel('tau');
% ylabel('b');

end

