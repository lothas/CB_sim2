%% Mixture of Experts
clear all; close all; clc;

%% Load data: 4 neuron CPG
load('MatsRandomRes_16_12_2016.mat','results','periods')
ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);
% For 4 neurons:
parametersCells = {'tau','b',...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};

targetCells = {'freq'};
[ sampl,targ ] = prepareData( results,periods,ids,parametersCells,targetCells,0 );

[sampl,targ ] = matsuoka_uniq(size(sampl,2),sampl,targ );

clear results periods ids_period ids_error ids
%% norm inputs - not necessary! NN toolbox already normalize!
% normParams = zeros(size(sampl, 1), 2);
% for i = 1:size(sampl, 1)
%     feat = sampl(i, :);
%     normParams(i, :) = [mean(feat), std(feat)];
%     sampl(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
%     
% end
% clear i feat normParams
%% Prepare the NN (="Experts")
learningRate = 0.5*10^-5;
decay = 0.9999;

expertCount = 5;      % how many "experts" (fitting NN)
numOfInputs = size(parametersCells,2); %how many inputs to each expert
maxEphocs = 200;      % max number of ephocs for each NN training
numOfIteretions = 50;  % number of loop interations

clusterCOG = cell(1,expertCount);  
dm_i = zeros(numOfInputs,expertCount);
expert_i_GroupSize = zeros(expertCount,numOfIteretions);
cluster_i__ind = cell(1,expertCount);
outMat = zeros(expertCount,size(sampl,2));
errMat = zeros(expertCount,size(sampl,2));
gateNN_perf_vec = zeros(1,numOfIteretions);
gateNet_targ = zeros(expertCount,size(sampl,2));
emptyGroupCounter = zeros(expertCount,1);

% initialize parameters, use first points for m
m = rand(expertCount, numOfInputs);

% define experts:
HiddenN = 10;
expertsNN = cell(2,expertCount);
for i=1:expertCount
    expertsNN{1,i} = feedforwardnet(HiddenN);
    expertsNN{1,i}.trainParam.showWindow = 0; % dont show training window
    expertsNN{1,i}.trainParam.epochs = maxEphocs;
end

% initial training (to initilazied the NN)
for j=1:expertCount
    initTrain = randsample(1:size(sampl,2),1000);
    [expertsNN{1,j}, expertsNN{2,j}] = train(expertsNN{1,j}, sampl(:,initTrain), targ(:,initTrain));
end

% define and initialize gate Network:
initTrain = randsample(1:size(sampl,2),10);
for j=1:expertCount
    tempNet = expertsNN{1,j};
    outMat(j,:) = tempNet(sampl);
    errMat(j,:) = outMat(j,:) - targ;
end
seMat = errMat.^2; % squar error
[~,best_expert_ind] = min(seMat,[],1);
for j=1:expertCount % find targets based on experts perf
    gateNet_targ(j,:) = (best_expert_ind == j);
end
% gate_init_targ = ones(expertCount,length(initTrain));
gateNet = patternnet(10); % view(gateNet);
gateNet.trainParam.showWindow = 0;
gateNet = train(gateNet,sampl,gateNet_targ);

%% Train the experts and the gate NN:
disp('runTime of training:');
tic
for i=1:numOfIteretions
    
    % devide to clusters with the gate network:
    gateOut = gateNet(sampl);
    gateOut_inx = vec2ind(gateOut); % indecies of points to #clusters
    
    % clustering to different experts
    for j=1:expertCount
        cluster_i__ind{1,j} = find(gateOut_inx == j);
        expert_i_GroupSize(j,i) = length(cluster_i__ind{1,j}); % check the size of each cluster
        
    end
    
    % training the experts;
    for j=1:expertCount
        tempNet = expertsNN{1,j};
        if expert_i_GroupSize(j,i) > 0 % only train if the cluster is not empty
            [expertsNN{1,j}, expertsNN{2,j}] = train(tempNet,...
                sampl(:,cluster_i__ind{1,j}), targ(:,cluster_i__ind{1,j}));
            % TODO: think what to do if the cluster is empty
        else
            % if the expert's cluster is empty, train the expert on the 'n'
            % points with the best probability
            [~,tempGroup] = sort(gateOut(j,:));
            best_tempGroup = tempGroup(1,1:1000);
            [expertsNN{1,j}, expertsNN{2,j}] = train(tempNet,...
                sampl(:,best_tempGroup), targ(:,best_tempGroup));
            emptyGroupCounter(j,1) = emptyGroupCounter(j,1) + 1;
        end
    end
    
    % run each expert on the entire data:
    for j=1:expertCount
        tempNet = expertsNN{1,j};
        outMat(j,:) = tempNet(sampl);
        errMat(j,:) = outMat(j,:) - targ;
    end
    seMat = errMat.^2; % squar error
    [~,best_expert_ind] = min(seMat,[],1);
    
    % retrain gate network:
    for j=1:expertCount % find targets based on experts perf
        gateNet_targ(j,:) = (best_expert_ind == j);
    end
    gateNet = train(gateNet,sampl,gateNet_targ);
    
end
toc

inx = vec2ind(gateNet(sampl));
    
for j=1:expertCount
   figure;
   groupInd = cluster_i__ind{1,j};
   tempNet = expertsNN{1,j};
   Outputs = tempNet(sampl);
   trOut = Outputs(:,groupInd);
   trTarg = targ(groupInd);
   plotregression(trTarg,trOut,'Train');
   clear groupInd tempNet Outputs trOut
end

expertsNames = cell(1,expertCount);
for j=1:expertCount
    expertsNames{1,j} = ['#',num2str(j),' expert'];
end
figure;
plot(1:numOfIteretions,expert_i_GroupSize); hold on;
xlabel('#iteretion');   ylabel('group size [#points]');
legend(expertsNames);

%% Visualized results of NN weights:
for j=1:expertCount
    NN_weights_matrix_plot(expertsNN{1,j},parametersCells);
end
%% clustering the data: devide for experts:
close all;
if false  % NN classifier for data gating: (pattern regongnition NN)
    percep_targ = zeros(expertCount,size(sampl,2));
    for j=1:expertCount
        percep_targ(j,:) = (best_expert_ind == j);
    end
    gateNet = patternnet(30); % view(gateNet);
    gateNet = train(gateNet,sampl,percep_targ);
    inx1 = gateNet(sampl);
%     [~,inx] = min(inx1,[],1);
    inx = vec2ind(inx1);
end

if false % NN classifier for data gating: (competetive NN)
    gateNet = competlayer(expertCount);
    gateNet = train(gateNet,sampl);
    inx1 = gateNet(sampl);
    inx = vec2ind(inx1);
end

if false % NN classifier for data gating gating using 'g':
    
    N = size(sampl,2);
    g = zeros(expertCount,N);
    for i=1:N
        input = sampl(:,i);
        g(:,i) = exp(m*input)/sum(exp(m*input));
%         g(:,i) = softmax(m*input);
        g(:,i) = round(g(:,i),4);
    end
    [~,inx] = max(g,[],1);
    clear N
end

bestExpertsInx = [inx;-best_expert_ind];
disp(['the number of train points correctly classify: ',...
    num2str(length(find(~(sum(bestExpertsInx,1))))),...
    ' out of ',num2str(length(bestExpertsInx))]);
%% check train results:
cluster_i_train_ind = cell(1,expertCount);
out_train = zeros(size(targ));
targ_train = [];
outM_train = [];
for j=1:expertCount
    cluster_i_train_ind{1,j} = find(inx == j);
    tempNet = expertsNN{1,j};
    out_train_temp = tempNet(sampl(:,cluster_i_train_ind{1,j}));
    targ_temp = targ(:,cluster_i_train_ind{1,j});
    targ_train = [targ_train,targ_temp];
    outM_train = [outM_train,out_train_temp];
    clear targ_temp out_train_temp
end

figure;
plotregression(targ_train,outM_train,'Train');

% calc R^2
err = targ_train-outM_train;
errVar = var(err,0,2);
inputVar = var(targ_train,0,2);

R_squar = 1-(errVar/inputVar);
disp(['the R^2 is: ',num2str(R_squar)]);

err = immse(outM_train,targ_train);
disp(['the MSE is: ',num2str(err)]);

%% for symmetric 2neuron CPG "knee checking"
% % run many genes with similar tau,T,b,c and different a and check the period
% % in this section I will check Matsuoka's approximation compare to the NN approximation.
% % given the same tau,b,c and different 'a'.
% 
% % attention change the genome file to match
% genome_file = 'MatsuokaGenome.mat';
% if false
%     nAnkle = 1;%1; % Number of ankle torques
%     nHip = 0;   % Number of hip torques
%     maxAnkle = 20;   % Max ankle torque
%     maxHip = 20;    % Max hip torque
%     Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
%     mamp = 0*Mamp;
%     N = nAnkle+nHip;
%     Mw = 10*ones(1,(2*N-1)*2*N);
%     mw = 0*Mw;
% % %    % for comparing a specific Matsuoka CPG to the NN
%     Keys = {'tau'   ,   'tav','tau_ratio', 'beta', 'amp_2n',    '2neuron_symm_weights' , 'ks_\tau',     'ks_c', 'IC_matsuoka';
%                1    ,       1,          1,      1,      2*N,                          1,        1 ,       2*N ,            0 };
%     Range = {  0.02 ,    0.01,        0.1,    0.2,     mamp,                        0.1,   -0.001 ,  -0.2*Mamp; % Min
%                    2,       2,         10,   10.0,     Mamp,                         10,   0.001 ,   0.2*Mamp}; % Max
% 
%     MutDelta0 = 0.04;   MutDelta1 = 0.02;
%     
%     save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
%         'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
%         'MutDelta0', 'MutDelta1', 'Keys', 'Range');
% end
%     
% close all; clc; %clear all;
% MML = MatsuokaML(); % calling Matsuoka class
% MML.perLim = [0.68 0.78];
% MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
% MML.tStep = 0.01;
% MML.tEnd = 30; % 15
% 
% N = 20;
% freq_Matsuoka = zeros(1,N);
% freq_fromNN = zeros(1,N);
% freq_fromCode = zeros(1,N);
% tau = 0.25;
% tau_ratio = 5;%2;
% T = tau*tau_ratio;
% b = 2.5;
% c = 5;
% a = linspace(1.6,3.4,N);
% tic
% for i=1:N
%     seq = [tau,T,tau_ratio,b,c,0,a(1,i),0,0,0];
%     [out, ~, ~] = MML.runSim(seq);
%         % Prepare output:
%     % Parameters
%     results(i).seq = seq;
%     results(i).tau = tau;
%     results(i).T = T;
%     results(i).b = b;
%     results(i).c = c;
%     results(i).a = a(1,i);
%     results(i).x0 = out.x0;
% 
%     % Results
%     results(i).periods = out.periods;
%     results(i).pos_work = out.pos_work;
%     results(i).neg_work = out.neg_work;
%     results(i).perError1 = out.perError1;
%     results(i).perOK1 = out.perOK1;
%     results(i).perError2 = out.perError2;
%     results(i).perOK2 = out.perOK2;
%     results(i).neuronActive = out.neuronActive;
%     results(i).neuronOsc = out.neuronOsc;
%     
%     freq_Matsuoka(1,i) = ((1/T)*sqrt(((tau+T)*b - a(1,i)*tau)/(a(1,i)*tau)))/(2*pi); % in [Hz]
%     
%     gateOut_2N = gateNet([tau;b;a(1,i);c]);
%     [~,maxGateOut_2N] = max(gateOut_2N,[],1);
%     tempNet = expertsNN{1,maxGateOut_2N};
%     freq_fromNN(1,i) = tempNet([tau;b;a(1,i);c]);
%     clear tempNet maxGateOut_2N gateOut_2N
%     
%     freq_fromCode(1,i)= 1/out.periods;
% end
% toc
% 
% figure;
% scatter(vertcat(results(:).a),freq_Matsuoka,'b','x'); hold on;
% scatter(vertcat(results(:).a),freq_fromNN,'r','d');
% scatter(vertcat(results(:).a),freq_fromCode,'g','o');
% xlabel('a');    ylabel('freq [Hz]'); grid on;
% legend('blue: Matsuoka estimation','red: freq from NN','green: freq from Code');
% title('frequency over a');