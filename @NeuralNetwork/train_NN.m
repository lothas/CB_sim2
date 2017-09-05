function obj = train_NN(obj)
% create and train the NN

% 'architecture' - row vector with the NN architecture
% 'inputs' - NN inputs
% 'targets' - NN targets

targets_num = size(obj.targets,1);

% Create and train the NN
% net = feedforwardnet(obj.architecture);
% net = feedforwardnet(obj.architecture,'useGPU','yes');
net = fitnet(obj.architecture);
% net = cascadeforwardnet(obj.architecture);

% net.trainFcn = obj.net_train_function;
% net.trainFcn = 'trainlm';
net.trainFcn = 'trainbr';
% net.trainFcn = 'trainscg';

% the two best options:
% net.trainFcn = 'trainrp';
% net.trainFcn = 'traingdm';

% dont show training window:
net.trainParam.showWindow = 1; 

if false % train only on some of the samples (not all of them)
    ind = randsample(size(obj.targets,2),50000);
    inputs = inputs(:,ind);
    targets = targets(:,ind);
end

% [net, tr] = train(net, obj.inputs, obj.targets,'useParallel','yes','useGPU','yes');
[net, tr] = train(net, obj.inputs, obj.targets);

% saving the 'net' object:
obj.net = net;
obj.net_train_perf = tr;

% calc NN output on the different groups
out_train = net(obj.inputs(:,tr.trainInd));
out_valid = net(obj.inputs(:,tr.valInd));
out_test = net(obj.inputs(:,tr.testInd));

% seperate the targets to the different groups
targ_train = obj.targets(:,tr.trainInd);
targ_valid = obj.targets(:,tr.valInd);
targ_test = obj.targets(:,tr.testInd);

% % % % test code, delete after use
% ind1 = find(obj.targets <= 2.5);
% out_test = net(obj.inputs(:,ind1));
% targ_test = obj.targets(:,ind1);
% % % % % % % % % % % % % % % % % % % % % 

% calculating NN perf for each target:
RMSE_train = zeros(1,targets_num);
RMSE_valid = zeros(1,targets_num);
RMSE_test = zeros(1,targets_num);

R2_train = zeros(1,targets_num);
R2_val = zeros(1,targets_num);
R2_test = zeros(1,targets_num);

slope_train = zeros(1,targets_num);
slope_val = zeros(1,targets_num);
slope_test = zeros(1,targets_num);
for i=1:targets_num
    [MSE_train,R2_train(1,i),slope_train(1,i)] =...
        obj.NN_perf_calc(targ_train(i,:),out_train(i,:),0);
    [MSE_val,R2_val(1,i),slope_val(1,i)] =...
        obj.NN_perf_calc(targ_valid(i,:),out_valid(i,:),0);
    [MSE_test,R2_test(1,i),slope_test(1,i)] =...
        obj.NN_perf_calc(targ_test(i,:),out_test(i,:),0);
    
    RMSE_train(1,i) = sqrt(MSE_train);
    RMSE_valid(1,i) = NaN;%sqrt(MSE_val);
    RMSE_test(1,i) = sqrt(MSE_test);
%     perf_train(1,i) = sqrt(immse(out_train(i,:),targ_train(i,:)));
%     perf_valid(1,i) = sqrt(immse(out_valid(i,:),targ_valid(i,:)));
%     perf_test(1,i) = sqrt(immse(out_valid(i,:),targ_valid(i,:)));
end

% save the perf in class:
obj.train_RMSE = RMSE_train;
obj.valid_RMSE = RMSE_valid;
obj.test_RMSE = RMSE_test;

obj.train_R2 = R2_train;
obj.valid_R2 = R2_val;
obj.test_R2 = R2_test;

obj.train_slope = slope_train;
obj.valid_slope = slope_val;
obj.test_slope = slope_test;

end



