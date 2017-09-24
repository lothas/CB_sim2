
%% NN testing Script:
% IMPORTANT: dont forget to load the right genome file and to uptade
% 'MatsuokaML.m' to the right rettings
% 

clear all; close all; clc;

%% Create genome (only if necessary)
genome_file = 'MatsuokaGenome_2Neuron_Symm.mat';
nAnkle = 1;%1; % Number of ankle torques
nHip = 0;   % Number of hip torques
maxAnkle = 10;   % Max ankle torque
maxHip = 10;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;

       % Narrow b narrow W narrow tau
Mw = 5;
mw = 0;
Keys = {'\tau_r', 'beta',     'amp_2n',    '2neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
              1 ,      1,            2,                         1,        1 ,          2,            0 };
Range = {  0.02 ,    0.2,        [0,0],                         0,   -0.001 ,[-0.2,-0.2]; % Min
           0.10  ,   2.5,      [10,10],                         5,    0.001 , [0.2,0.2]}; % Max

MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

clear all
%%

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau' ,'b', 'c', 'NR', 'a',...
    'k_tau','k_{c}'};
% "NR" - not relevnt param 

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;

% % file name for uploading:
% results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_only_osc_1.mat'};

results_fileName = {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_only_osc_1.mat',...
    'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_only_osc_2.mat'};

%% Load data:
results = load_results(results_fileName);

%% plot example:
clc; close all

N = length(results);
rand_id = randsample(1:N,1);

[out, ~, signal] = MML.runSim(results(rand_id).seq);

figure;
subplot(2,1,1);
plot(signal.T,signal.X);
xlabel('time[sec]');    ylabel('X_i');
title({'X_i over time',...
    ['id #',num2str(rand_id),...
    '    periods: ',...
    num2str(results(rand_id).periods(1))]});
subplot(2,1,2)
plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');

clear out signal N rand_id

%% get and filter periods:
% % get oscillating:
% [results,periods,seq,~] = get_CPGs(results,'osc','2N',MML);

% % % get oscillating in period range:
[results,periods,seq,~] = get_CPGs(results,'osc_in_per_range','2N',MML);
% 
figure;
boxplot(seq','orientation','horizontal','labels',seqOrder)

plot_param_hist(seq,periods,seqOrder)

%% remove outliers:
% https://www.mathworks.com/matlabcentral/answers/121247-how-can-i-detect-and-remove-outliers-from-a-large-dataset

ids = true(1,size(seq,2));

for i=1:(MML.Gen.Length+1)
    
    if i > MML.Gen.Length
        vector = periods;
    else
        vector = seq(i,:);
    end

    percntiles = prctile(vector,[5 95]); %5th and 95th percentile
    % the distance between the %5th and 95th percentiles is four stdevs
    
    outlierIndexes = vector < percntiles(1) | vector > percntiles(2);
    
    ids = ~outlierIndexes & ids;
%     % Extract outlier values:
%     outliers = vector(outlierIndexes);
%     % Extract non-outlier values:
%     nonOutliers = vector(~outlierIndexes);
end

periods = periods(:,ids);
seq = seq(:,ids);

figure;
subplot(2,1,1);
boxplot(seq','orientation','horizontal','labels',seqOrder);
subplot(2,1,2);
boxplot(periods,'orientation','horizontal','labels',{'periods'});

clear ids
%% Prepare NN inputs and outputs:
% input_names = {'b','tau','a'};
% output_names = {'periods'};

input_names = {'tau','a','periods'};
output_names = {'b'};

[sampl,targ] = ...
    prepare_NN_data(input_names,output_names,...
    seqOrder,seq,periods);

%% some 2D histograms:
close all;

figure;
subplot(2,2,1);
histogram2(sampl(1,:),targ(1,:),100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of \tau and b');
xlabel('tau');
ylabel('b');

subplot(2,2,2);
histogram2(sampl(3,:),targ(1,:),100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of period and b');
xlabel('period');
ylabel('b');

subplot(2,2,3);
histogram2(sampl(2,:),targ(1,:),100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title({'2D distribution of \tau and a',...
    'NOTE: 1+tau/T < a < 1+b'});
xlabel('a');
ylabel('b');

figure;
subplot(2,2,1);
histogram2(sampl(3,:),sampl(1,:),100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of \tau and the periods');
xlabel('period');
ylabel('tau');

subplot(2,2,2);
histogram2(sampl(2,:),sampl(1,:),100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of \tau and "a"');
xlabel('a');
ylabel('tau');

subplot(2,2,3);
histogram2(targ(1,:),sampl(1,:),100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of \tau and "b"');
xlabel('b');
ylabel('tau');

figure;
subplot(2,2,1);
histogram2(sampl(2,:),sampl(1,:),100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of "a" and \tau');
xlabel('a');
ylabel('tau');

subplot(2,2,2);
histogram2(sampl(2,:),sampl(3,:),100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of "a" and the periods');
xlabel('a');
ylabel('periods');

subplot(2,2,3);
histogram2(sampl(2,:),targ(1,:),100,...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('2D distribution of "a" and "b"');
xlabel('a');
ylabel('b');

%%
close all;

Title = '3D histogram '; 
XYbinNum = 20;
zbinNum = 3;

Labels = {'\tau','a','periods'};
hist3D_CS(sampl(1,:),sampl(2,:),sampl(3,:),Labels,Title,XYbinNum,zbinNum);

Labels = {'b','a','periods'};
hist3D_CS(targ(1,:),sampl(2,:),sampl(3,:),Labels,Title,XYbinNum,zbinNum);

Labels = {'b','a','tau'};
hist3D_CS(targ(1,:),sampl(2,:),sampl(1,:),Labels,Title,XYbinNum,zbinNum);
%% Neural Network #1a: (CPU training)
architecture = [30];

net = fitnet(architecture);
% net = feedforwardnet(architecture);
net.trainFcn = 'trainbr';
net.trainParam.showWindow = 1; 
net.trainParam.showCommandLine = 1;
net.trainParam.epochs = 1000;
% net.performParam.normalization = 'percent'; %It must be 'none', 'standard' or 'percent'

[net, tr] = train(net, sampl, targ);
[r,m,b] = regression(targ,net(sampl));
disp(['r = ',num2str(r),'    m = ',num2str(m),'    b = ',num2str(b)]);
clear r m b

figure;
histogram(targ,100,'Normalization','pdf'); hold on;
histogram(net(sampl),100,'Normalization','pdf');
legend('targets','NN outputs');

% [net, tr] = train(net, sampl_n, targ_n);
% figure;
% histogram(targ_n,100,'Normalization','pdf'); hold on;
% histogram(net(sampl_n),100,'Normalization','pdf');
% legend('targets','NN outputs');
%% Neural Network #1b: (GPU training)
architecture = [10,20,20,20,5];

net = fitnet(architecture);
% net = feedforwardnet(architecture);
net.trainParam.showWindow = 1; 
net.trainParam.showCommandLine = 1;
net.trainParam.epochs = 1000;
% net.performParam.normalization = 'percent'; %It must be 'none', 'standard' or 'percent'

[net, tr] = train(net, sampl, targ,'useGPU','yes');
[r,m,b] = regression(targ,net(sampl));
disp(['r = ',num2str(r),'    m = ',num2str(m),'    b = ',num2str(b)]);
clear r m b

figure;
histogram(targ,100,'Normalization','pdf'); hold on;
histogram(net(sampl),100,'Normalization','pdf');
legend('targets','NN outputs');

%%
clc;% close all
% 
% input_names = {'tau','a','periods'};
% output_names = {'b'};
% plot_param_hist(seq,periods,seqOrder)

N = 100;
M = 6;

tau = 0.03;
a = linspace(0.1,4.9,M);
p = linspace(0.1,1,N);%0.7;

NNout = zeros(M,N);

for j=1:M
    for i=1:N
        NNin = [tau;a(1,j);p(1,i)];
        NNout(j,i) = net(NNin);
%         NNin_n = mapstd('apply',NNin,sampl_s);
%         NNout_n = sim(net,NNin_n);
%         NNout(j,i) = mapstd('reverse',NNout_n,targ_s);        
    end
end

figure; hold on;
for j=1:M
    plot(p,NNout(j,:),'Marker','o');
    LABELS{1,j} = sprintf('a = %f',a(1,j));
end
grid minor;
legend(LABELS);
title(['tau = ',num2str(tau)]);
xlabel('periods');
ylabel('NN outputs');

clear N tau a p NNout LABELS rand_id N
%% % test NN with training data (only with des period)
close all

desPeriod = MML.perLim(1) + ...
                 rand(1,size(targ,2))*(MML.perLim(2)-MML.perLim(1));
% % desPeriod = 1 + rand(1,size(targ,2))*(5-1);
% 
[NN_in,~] = ...
    prepare_NN_data(input_names,output_names,...
    seqOrder,seq,desPeriod);

% % % %  test on rande seq #1:
% rand_seq = MML.Gen.RandSeq(length(targ));
% [NN_in,~] = ...
%     prepare_NN_data(input_names,output_names,...
%     seqOrder,rand_seq',desPeriod);

% % % get NN output:
NN_out = net(NN_in);
plotregression(targ,NN_out)

% % % % % find nearest neighbours:
[IDX,D] = knnsearch(NN_out',targ');

figure;
scatter(NN_out(1:1000),targ(1,IDX(1:1000)));
xlabel('NN output "b"');
ylabel('nearest neighbor from the training data')
% axis([0 2.5 0 2.5])
% % % % % % % % % % % % % % % 

% NN_out = net(sampl(:,tr.testInd));
% plotregression(targ(:,tr.testInd),NN_out)

figure;
histogram(NN_out,100,'Normalization','pdf');
title('histogram of NN_{output}');
xlabel('NN_{output}');

figure;
histogram(periods,100,'Normalization','pdf');
title('histogram of periods');
xlabel('periods');

figure;
subplot(2,1,1);
boxplot(sampl','orientation','horizontal','labels',input_names)
subplot(2,1,2);
boxplot(targ','orientation','horizontal','labels',output_names)

clear NN_in NN_out desPeriod