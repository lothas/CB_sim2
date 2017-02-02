%% Mixture of Experts
clear all; close all; clc;

%% Load data: 4 neuron CPG
howMuchData = 600000;
disp(['Loading approx ',num2str(howMuchData),' samples...']);
[results,periods]=load_data_4N_CPG(howMuchData);

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

%% Load Data: 2neurons symmetric case
howMuchData = 400000;
[results,periods]=load_data_2N_Symm_CPG(howMuchData);

% identify which neurons have periods
ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);

% prepare NN inputs and outputs:
parametersCells = {'tau','b','a','s'};
targetCells = {'freq'};
[ sampl,targ ] = prepareData_2Neurons( results,periods,ids,parametersCells,targetCells);

clear ids_period ids_error

%% Load Data: 2neurons general case
howMuchData = 400000;
[results,periods]=load_data_2N_Symm_CPG(howMuchData);

% identify which neurons have periods
ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);

parametersCells = {'tau','b','w12','w21'};
targetCells = {'freq'};
sampl = vertcat(results(ids).seq)'; % extracting the parameters from the structure to a matrix
sampl = sampl([1,2,5,6],:); % get rid of the the unimportant 'gen' (the k's for the control).
targ= 1./periods(ids);

clear ids_period ids_error nSims

%% norm inputs - not necessary! NN toolbox already normalize!
normParams = zeros(size(sampl, 1), 2);
for i = 1:size(sampl, 1)
    feat = sampl(i, :);
    normParams(i, :) = [mean(feat), std(feat)];
    sampl(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
    
end
clear i feat normParams
%%
expertCount = 3;      % how many "experts" (fitting NN)
numOfInputs = size(parametersCells,2); %how many inputs to each expert
maxEphocs = 5;      % max number of ephocs for each NN training
numOfIteretions = 10;  % number of loop interations
ExpertHidLayer = 1; % num of hidden layer in each expert
ExpertHidNueron = 10; % num of neurons in each hidden layer
GateHidLayer = 1; % num of hidden layer in gateNN
GateHidNueron = 10; % num of neurons in each hidden layer
competetiveFlag = 3; % if '1'- "winner takes all"
                     %    '2'- "chance for everybody"
                     %    '3'- out = expertsOut * gateOut
[ expertsNN,gateNet,expert_i_GroupSize,gateNN_perf_vec,Experts_perf_mat,Moe_perf_over_iter,emptyGroupIndecator ] = ...
    my_MoE_train(sampl,targ,expertCount,numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
                GateHidLayer,GateHidNueron,competetiveFlag);

[netOut,gateOut,targ,~,cluster_i_train_ind] = my_MoE_testNet(sampl,targ,expertsNN,...
    gateNet,competetiveFlag);

my_MoE_plotPerf(netOut,targ,gateOut,cluster_i_train_ind,Moe_perf_over_iter,...
    gateNN_perf_vec,expert_i_GroupSize,Experts_perf_mat,emptyGroupIndecator,...
    'both',competetiveFlag);


[~,~] = NN_perf_calc(targ,netOut,1,0);

%% test group analysis:
load('MatsRandomRes_21_12_2016.mat','results','periods');
ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);
[ test_sampl,test_targ ] = prepareData( results,periods,ids,parametersCells,targetCells,0 );
[test_sampl,test_targ ] = matsuoka_uniq(size(test_sampl,2),test_sampl,test_targ );
clear results periods ids_period ids_error ids

[netOut_test,gateOut_test,targ_test,~,cluster_i_train_ind] = my_MoE_testNet(test_sampl,test_targ,expertsNN,...
    gateNet,competetiveFlag);

my_MoE_plotPerf(netOut_test,targ_test,gateOut,cluster_i_train_ind,competetiveFlag);

% [errMSE,R_squar] = NN_perf_calc(targ_test,netOut_test,dispFlag,GraphicsFlag);
%% Visualized results of NN weights:
for j=1:expertCount
    NN_weights_matrix_plot(expertsNN{1,j},parametersCells);
end

%% many training cross validation
% compare two methodes of my MoE code. 1st method is completely competetive
% (winner takes all) and th 2nd one is more probabilistic.
numOfrepeats = 5;
numOfIteretions = 50;  % number of loop interations
expertCount = [2,3,4,5];      % how many "experts" (fitting NN)
competetiveFlag_vec = [1,2,3];
R_squar_mat = zeros(length(competetiveFlag_vec),length(expertCount),numOfrepeats);
errMSE_mat = zeros(length(competetiveFlag_vec),length(expertCount),numOfrepeats);

for k=1:length(competetiveFlag_vec)
    competetiveFlag = competetiveFlag_vec(1,k);
    for n=1:length(expertCount)
        for i=1:numOfrepeats
            [expertsNN,gateNet,~,~,~,~,~] = my_MoE_train(sampl,targ,expertCount(1,n),...
                numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
                        GateHidLayer,GateHidNueron,competetiveFlag);
                    
            [netOutt,gateOut,targ_temp,~,~] = my_MoE_testNet(sampl,targ,expertsNN,gateNet,competetiveFlag);
            
            [errMSE_mat(k,n,i),R_squar_mat(k,n,i)] = NN_perf_calc(targ_temp,netOut,0,0);                    
        end
    end
end

% plotting bar graphs of avarge MSE + errorbars
meanMse_CF1 = mean(R_squar_mat(1,:,:),2); % "CF"=Competetive Flag
stdMse_CF1 = std(errMSE_mat(1,:,:),0,2);
meanMse_CF2 = mean(R_squar_mat(1,:,:),2); 
stdMse_CF2 = std(errMSE_mat(1,:,:),0,2);
meanMse_CF3 = mean(R_squar_mat(1,:,:),2);
stdMse_CF3 = std(errMSE_mat(1,:,:),0,2);

stdevs = [stdMse_CF1;stdMse_CF2;stdMse_CF3];
means = [meanMse_CF1;meanMse_CF2;meanMse_CF3];

Names={' ';'2 Experts';' ';'3 Experts';'  ';'4 Experts';' ';'5 Experts'};
label_Y = {'MSE'};
graph_title = {'performance over different #experts and different methods'};
graph_legend = {'winner takes all','chance for everybody','out = expertsOut * gateOut'};

plotBars_with_errBars( means',stdevs',Names,label_Y,graph_title,graph_legend)
%% regression graphs for the paper's algorithm:
expertCount = 10;
numOfIteretions = 1000;

[ExpertsWeights, gateWeights,errs] = paper_MoE_train(sampl, targ, expertCount, numOfIteretions,0.001, 0.995);

[~, freq_fromNN,g] = paper_MoE_test(sampl,targ, ExpertsWeights, gateWeights,1);

paper_MoE_plotPerf(sampl,targ,freq_fromNN,g,errs);

%% for symmetric 2neuron CPG "knee checking"
% run many genes with similar tau,T,b,c and different a and check the period
% in this section I will check Matsuoka's approximation compare to the NN approximation.
% given the same tau,b,c and different 'a'.
% attention change the genome file to match!!!

close all; clc; %clear all;
% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01;
MML.tEnd = 30; % 15

N = 30;

% CPG parameters:
tau = 0.25;
tau_ratio = 5;%2;
T = tau*tau_ratio;
b = 2.5;
c = 5;
a = linspace(1.6,3.4,N);

% MoE parameters:
paper_MoE_flag = true;
competetiveFlag = 1;
expertCount = 3;%[2,3,4]; 
numOfIteretions = 10;

if false
    genome_file = 'MatsuokaGenome.mat';
    nAnkle = 1;%1; % Number of ankle torques
    nHip = 0;   % Number of hip torques
    maxAnkle = 20;   % Max ankle torque
    maxHip = 20;    % Max hip torque
    Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
    mamp = 0*Mamp;
    N = nAnkle+nHip;
    Mw = 10*ones(1,(2*N-1)*2*N);
    mw = 0*Mw;
% %    % for comparing a specific Matsuoka CPG to the NN
    Keys = {'tau'   ,   'tav','tau_ratio', 'beta', 'amp_2n',    '2neuron_symm_weights' , 'ks_\tau',     'ks_c', 'IC_matsuoka';
               1    ,       1,          1,      1,      2*N,                          1,        1 ,       2*N ,            0 };
    Range = {  0.02 ,    0.01,        0.1,    0.2,     mamp,                        0.1,   -0.001 ,  -0.2*Mamp; % Min
                   2,       2,         10,   10.0,     Mamp,                         10,   0.001 ,   0.2*Mamp}; % Max

    MutDelta0 = 0.04;   MutDelta1 = 0.02;
    
    save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
        'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
        'MutDelta0', 'MutDelta1', 'Keys', 'Range');
end   % define the genome file structure for matsuokaML.runSim

disp('start with the sim:');
freq_Matsuoka = zeros(1,N);
freq_fromCode = zeros(1,N);
for i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
    disp(['at sim #',num2str(i)]);
    seq = [tau,T,tau_ratio,b,c,0,a(1,i),0,0,0];
    [out, ~, ~] = MML.runSim(seq);
        % Prepare output:
    % Parameters
    results_2N_sim(i).seq = seq;
    results_2N_sim(i).tau = tau;
    results_2N_sim(i).T = T;
    results_2N_sim(i).b = b;
    results_2N_sim(i).c = c;
    results_2N_sim(i).a = a(1,i);
    results_2N_sim(i).x0 = out.x0;

    % Results
    results_2N_sim(i).periods = out.periods;
    results_2N_sim(i).pos_work = out.pos_work;
    results_2N_sim(i).neg_work = out.neg_work;
    results_2N_sim(i).perError1 = out.perError1;
    results_2N_sim(i).perOK1 = out.perOK1;
    results_2N_sim(i).perError2 = out.perError2;
    results_2N_sim(i).perOK2 = out.perOK2;
    results_2N_sim(i).neuronActive = out.neuronActive;
    results_2N_sim(i).neuronOsc = out.neuronOsc;

    freq_Matsuoka(1,i) = ((1/T)*sqrt(((tau+T)*b - a(1,i)*tau)/(a(1,i)*tau)))/(2*pi); % in [Hz]

    freq_fromCode(1,i)= 1/out.periods;
end 
disp('sim end...');
%% for symmetric 2neuron CPG "knee checking"
% check MoE results, 'tau','b' and 'c' fixed 
% 'a' is changing within a certain range
expertsWeights_store = cell(1,length(expertCount));
gateWeights_store = cell(1,length(expertCount));
g = cell(1,length(expertCount));
freq_fromNN = zeros(length(expertCount),N);
belongToExpert_vec = zeros(length(expertCount),N);

for n=1:length(expertCount)   
    switch paper_MoE_flag
        case 1 %train the MoE with the paper's method
            [ExpertsWeights, gateWeights] = paper_MoE_train(sampl, targ, expertCount(1,n), numOfIteretions, 0.001, 0.9995);
            expertsWeights_store{1,n} = ExpertsWeights;
            gateWeights_store{1,n} = gateWeights;
                    v = [ones(1,N)*tau; ones(1,N)*b; a; ones(1,N)*c];
            [~, freq_fromNN(n,:),~] = paper_MoE_test(v, (1./horzcat(results_2N_sim(:).periods)), ExpertsWeights, gateWeights,0);
            clear v
        case 0  %train the MoE with our method
            [ expertsNN,gateNet,~,~,~,Moe_perf_over_iter,~]= ...
                my_MoE_train(sampl,targ,expertCount,numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
                             GateHidLayer,GateHidNueron,competetiveFlag);
            for i=1:N
                targ_temp = (1./horzcat(results_2N_sim(:).periods));
                [freq_fromNN(n,i),~,~,~,~] =...
                    my_MoE_testNet(sampl,targ_temp,expertsNN,gateNet,competetiveFlag);
                clear targ_temp        
            end
    end
end

figure;
plot(vertcat(results_2N_sim(:).a),freq_fromCode,'g-o'); hold on
plot(vertcat(results_2N_sim(:).a),freq_Matsuoka,'b-o');
plot(vertcat(results_2N_sim(:).a),freq_fromNN);
xlabel('a');    ylabel('freq [Hz]'); grid on;
legendNames = {'real freq','Matsuoka estimation',cell(1,length(expertCount))};
for j=3:(length(expertCount)+2)
    legendNames{1,j} = ['#',num2str(expertCount(1,j-2)),' expert'];
end
legend(legendNames);
title('frequency over "a" with different #experts');

[~,~] = NN_perf_calc(freq_fromCode,(freq_fromNN(n,:)),1,1);

legendNames = cell(1,expertCount);
for j=1:expertCount
    legendNames{1,j} = ['#',num2str(j),' expert'];
end
    
figure;
if paper_MoE_flag
    subplot(2,1,1);
    plot(vertcat(results_2N_sim(:).a),freq_fromCode,'g-o'); hold on
    plot(vertcat(results_2N_sim(:).a),freq_Matsuoka,'b-o');
    plot(vertcat(results_2N_sim(:).a),freq_fromNN(n,:),'r-o'); grid minor;
    xlabel('a'); ylabel('\omega_{n}');
    title('frequency over "a"');

    subplot(2,1,2);
    bar(g','stacked'); xlabel('#sample');
    legend(legendNames);
else
    plot(vertcat(results_2N_sim(:).a),freq_fromCode,'g-o'); hold on
    colors = rand(expertCount,3);
    for j=1:N
        expertNum = belongToExpert_vec(n,j);
        h = plot(vertcat(results_2N_sim(j).a),freq_fromNN(n,j),'Color',colors(expertNum,:),'LineStyle','none');
        h.Marker = 'o';
    end
    xlabel('a'); ylabel('\omega_{n}');
    legend(legendNames)
    title('frequency over "a"');
end


