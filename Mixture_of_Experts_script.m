%% Mixture of Experts
clear all; close all; clc;

%% Load data: 4 neuron CPG
load('MatsRandomRes_16_12_2016.mat','results','periods')
results1 = results; periods1 = periods;
clear results periods
load('MatsRandomRes_18_12_2016.mat','results','periods')
results2 = results; periods2 = periods;
clear results periods
load('MatsRandomRes_19_12_2016.mat','results','periods')
results3 = results; periods3 = periods;
clear results periods
load('MatsRandomRes_20_12_2016.mat','results','periods')
results4 = results; periods4 = periods;
clear results periods
results = horzcat(results1,results2,results3,results4);
periods = horzcat(periods1,periods2,periods3,periods4);
clear results1 results2 results3 results4...
    periods1 periods2 periods3 periods4

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
maxEphocs = 200;      % max number of ephocs for each NN training
numOfIteretions = 50;  % number of loop interations
ExpertHidLayer = 1; % num of hidden layer in each expert
ExpertHidNueron = 5; % num of neurons in each hidden layer
GateHidLayer = 1; % num of hidden layer in gateNN
GateHidNueron = 5; % num of neurons in each hidden layer
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
[test_sampl,test_targ ] = matsuoka_uniq(size(test_sampl,2),test_sampl,test_targ );
clear results periods ids_period ids_error ids

[R_squar_test,errMSE_test] = my_MoE_testing(test_sampl,test_targ,expertsNN,...
    gateNet,GraphicsFlag,competetiveFlag);

%% many training cross validation
% compare two methodes of my MoE code. 1st method is completely competetive
% (winner takes all) and th 2nd one is more probabilistic.
numOfrepeats = 5;
GraphicsFlag = 0;
numOfIteretions = 50;  % number of loop interations
expertCount = [2,3,4,5];      % how many "experts" (fitting NN)
R_squar_vec_competetive = zeros(numOfrepeats,1);
errMSE_vec_competetive = zeros(numOfrepeats,1);
R_squar_vec_prob = zeros(numOfrepeats,length(expertCount));
errMSE_vec_prob = zeros(numOfrepeats,length(expertCount));
competetiveFlag = 1;
for n=1:length(expertCount)
    for i=1:numOfrepeats
        [~,~,~,~,~,~,R_squar_vec_competetive(i,n),errMSE_vec_competetive(i,n),~]...
            = my_MoE_train(sampl,targ,expertCount(1,n),...
        numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
        GateHidLayer,GateHidNueron,GraphicsFlag,competetiveFlag);
    end

    competetiveFlag = 0; 
    for i=1:numOfrepeats
        [~,~,~,~,~,~,R_squar_vec_prob(i,n),errMSE_vec_prob(i,n),~]...
            = my_MoE_train(sampl,targ,expertCount(1,n),...
        numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
        GateHidLayer,GateHidNueron,GraphicsFlag,competetiveFlag);
    end
end
% plotting bar graphs of avarge MSE + errorbars
meanMse_comp = mean(errMSE_vec_competetive,1);
stdMse_comp = std(errMSE_vec_competetive,0,1);
meanMse_prob = mean(errMSE_vec_prob,1);
stdMse_prob = std(errMSE_vec_prob,0,1);
stdevs = [stdMse_comp;stdMse_prob];
means = [meanMse_comp;meanMse_prob];
% means = means';     stdevs = stdevs';
Names={' ';'2 Experts';' ';'3 Experts';'  ';'4 Experts';' ';'5 Experts'};
numgroups = size(means, 2); 
numbars = size(means, 1); 
groupwidth = min(0.8, numbars/(numbars+1.5));
figure; hold on
h=bar(means');
set(gca,'XTickLabel',Names);
set(h,'BarWidth',1);
for k=1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*k-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x, means(k,:), stdevs(k,:), 'k', 'linestyle', 'none');
end
legend('competetive','prob');
ylabel('MSE');
title('MoE with different amount of Experts');
hold off

%% regression graphs for the paper's algorithm:
expertCount = 100;
numOfIteretions = 5;

[ExpertsWeights, gateWeights] = paper_MoE_train(sampl, targ, expertCount, numOfIteretions, 0.001, 0.95);

[~, freq_fromNN] = paper_MoE_test(sampl,targ, ExpertsWeights, gateWeights);

% [ExpertsWeights, gateWeights] = paper_MoE_train_nonlinear(sampl, targ, expertCount, numOfIteretions, 0.0001, 0.9995);
% 
% [~, freq_fromNN] = paper_MoE_test_nonlinear(sampl,targ, ExpertsWeights, gateWeights);

[errMSE,R_squar] = NN_perf_calc(targ,freq_fromNN',1,1);

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

N = 20;

% CPG parameters:
tau = 0.25;
tau_ratio = 5;%2;
T = tau*tau_ratio;
b = 2.5;
c = 5;
a = linspace(1.6,3.4,N);

% MoE parameters:
paper_MoE_flag = true;
expertCount = 20;%[2,3,4,5]; 
numOfIteretions = 100;

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

for n=1:length(expertCount)
    disp('start with the sim:');
    freq_Matsuoka = zeros(1,N);
    freq_fromCode = zeros(1,N);
    for i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation) 
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
      
    if paper_MoE_flag %train the MoE and
        [ExpertsWeights, gateWeights] = paper_MoE_train(sampl, targ, expertCount(1,n), numOfIteretions, 0.00001, 0.9995);
    else
        [ expertsNN,gateNet,~,~,~,~,~,~,~] = my_MoE_train(sampl,targ,...
        expertCount(1,n),numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
        GateHidLayer,GateHidNueron,GraphicsFlag,competetiveFlag);
    end
    %
    freq_fromNN = zeros(1,N);
    if paper_MoE_flag
        v = [ones(1,N)*tau;
            ones(1,N)*b;
            a;
            ones(1,N)*c];
        [~, freq_fromNN] = paper_MoE_test(v, (1./horzcat(results_2N_sim(:).periods)), ExpertsWeights, gateWeights);
        clear v
    else
        for i=1:N
                gateOut_2N = gateNet([tau;b;a(1,i);c]);
                [~,maxGateOut_2N] = max(gateOut_2N,[],1);
                tempNet = expertsNN{1,maxGateOut_2N};
                freq_fromNN(1,i) = tempNet([tau;b;a(1,i);c]);
    %             plotMarkets = {'b-o','r-o','g-o'};
    %             h = plot(a(1,i),freq_fromNN(1,i),plotMarkets(1,maxGateOut_2N));
            clear tempNet maxGateOut_2N gateOut_2N           
        end
    end
    
    figure;
    plot(vertcat(results_2N_sim(:).a),freq_Matsuoka,'b-o'); hold on;
    plot(vertcat(results_2N_sim(:).a),freq_fromNN,'r-o');
    plot(vertcat(results_2N_sim(:).a),freq_fromCode,'g-o');
    xlabel('a');    ylabel('freq [Hz]'); grid on;
    legend('Matsuoka estimation','freq from NN','freq from Code');
    title('frequency over a');
    
end