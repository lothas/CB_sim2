%% Mixture of Experts
clear all; close all; clc;

%% Load data: 4 neuron CPG
load('MatsRandomRes_16_12_2016.mat','results','periods')
results1 = results; periods1 = periods;
clear results periods
load('MatsRandomRes_18_12_2016.mat','results','periods')
results2 = results; periods2 = periods;
clear results periods
results = horzcat(results1,results2);
periods = horzcat(periods1,periods2);
clear results1 results2 periods1 periods2

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
%%
expertCount = 3;      % how many "experts" (fitting NN)
numOfInputs = size(parametersCells,2); %how many inputs to each expert
maxEphocs = 200;      % max number of ephocs for each NN training
numOfIteretions = 50;  % number of loop interations
ExpertHidLayer = 1; % num of hidden layer in each expert
ExpertHidNueron = 10; % num of neurons in each hidden layer
GateHidLayer = 1; % num of hidden layer in gateNN
GateHidNueron = 10; % num of neurons in each hidden layer
GraphicsFlag = 1; % if to plot or not
competetiveFlag = 0; %if 1'- than the clustering is being done
%                       by highest P takes all. if not '0' than
%                       each sample as a chance to go to each cluster based
%                       on the probability from the gate network
[ expertsNN,gateNet,gateNet_perf,...
    expert_i_GroupSize,gateNN_perf_vec,Experts_perf_mat,...
    R_squar,errMSE,emptyGroupIndecator ] = ...
    my_MoE_train(sampl,targ,expertCount,...
    numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
    GateHidLayer,GateHidNueron,GraphicsFlag,competetiveFlag);

%% Visualized results of NN weights:
for j=1:expertCount
    NN_weights_matrix_plot(expertsNN{1,j},parametersCells);
end

%% more plots:
expertsNames = cell(1,expertCount);
for j=1:expertCount
    expertsNames{1,j} = ['#',num2str(j),' expert'];
end
figure;
subplot(2,2,1)
plot(1:numOfIteretions,expert_i_GroupSize,'-o'); hold on;
title('cluster size over #iteration');
xlabel('#iteretion');   ylabel('group size [#points]');
legend(expertsNames);

subplot(2,2,2) % expert perf over #interation
plot(1:numOfIteretions,Experts_perf_mat,'-o'); hold on;
title('expert MSE over #interation');
xlabel('#iteretion');   ylabel('performance [MSE]');

subplot(2,2,3) % expert perf over #interation
plot(1:numOfIteretions,double(emptyGroupIndecator),'-o'); hold on;
title('indication on empty clusters: "1" means empty');
xlabel('#iteretion');   ylabel('"1"=empty   "0"-not empty');

subplot(2,2,4) % gateNet perf over #interation
plot(1:numOfIteretions,gateNN_perf_vec,'-o'); hold on;
title('gateNet perf (crossEntropy) over #interation');
xlabel('#iteretion');   ylabel('performance [MSE]');

%% test group analysis:
load('MatsRandomRes_21_12_2016.mat','results','periods');
ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);
[ test_sampl,test_targ ] = prepareData( results,periods,ids,parametersCells,targetCells,0 );
[test_sampl,test_targ ] = matsuoka_uniq(size(test_sampl,2),test_sampl,targ );
clear results periods ids_period ids_error ids

[R_squar_test,errMSE_test] = my_MoE_testing(test_sampl,test_targ,expertsNN,...
    gateNet,GraphicsFlag,competetiveFlag);

%% many training cross validation
GraphicsFlag = 0;
for i=1:10
    [~,~,~,~,~,~,R_squar,errMSE,~] = my_MoE_train(sampl,targ,expertCount,...
    numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
    GateHidLayer,GateHidNueron,GraphicsFlag,competetiveFlag);
end
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