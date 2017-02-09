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
howMuchData = 10000;
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
%% our MoE
expertCount = 2;      % how many "experts" (fitting NN)
numOfInputs = size(parametersCells,2); %how many inputs to each expert
maxEphocs = 5;      % max number of ephocs for each NN training
numOfIteretions = 100;  % number of loop interations
ExpertHidLayer = 1; % num of hidden layer in each expert
ExpertHidNueron = 10; % num of neurons in each hidden layer
GateHidLayer = 1; % num of hidden layer in gateNN
GateHidNueron = 10; % num of neurons in each hidden layer
competetiveFlag = 3; % if '1'- "winner takes all"
                     %    '2'- "chance for everybody"
                     %    '3'- out = expertsOut * gateOut
[ expertsNN,gateNet,expert_i_GroupSize,gateNN_perf_vec,Experts_perf_mat,Moe_perf_over_iter,emptyGroupIndecator,~,~ ] = ...
    my_MoE_train(sampl,targ,expertCount,numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
                GateHidLayer,GateHidNueron,competetiveFlag);

[netOut,gateOut,targ,~,cluster_i_train_ind] = my_MoE_testNet(sampl,targ,expertsNN,...
    gateNet,competetiveFlag);

my_MoE_plotPerf(netOut,targ,gateOut,cluster_i_train_ind,Moe_perf_over_iter,...
    gateNN_perf_vec,expert_i_GroupSize,Experts_perf_mat,emptyGroupIndecator,...
    'both',competetiveFlag);


[~,~] = NN_perf_calc(targ,netOut,1,0,'train');

%% using the paper script
expertCount = 2;
numOfIteretions = 50;

[ExpertsWeights, gateWeights,errs,~,~] = paper_MoE_train(sampl, targ, expertCount, numOfIteretions,0.001, 0.995);

[~, freq_fromNN,g] = paper_MoE_test(sampl,targ, ExpertsWeights, gateWeights,1);

paper_MoE_plotPerf(sampl,targ,freq_fromNN,g,errs);

%% Visualized results of NN weights:
for j=1:expertCount
    NN_weights_matrix_plot(expertsNN{1,j},parametersCells);
end

%% define Matsuoka Class and the "Genome file":
close all; clc; 

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
% %     % 2neuron symmetric specific range%%
Keys = {'\tau_r', 'beta',     'amp_2n',        '2neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
              1 ,      1,          2*N,                             1,        1 ,       2*N ,            0 };
Range = {  0.02 ,      1,         mamp,                             1,   -0.001 ,  -0.2*Mamp; % Min
           0.6  ,    8.0,         Mamp,                             6,   0.001 ,   0.2*Mamp}; % Max

MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');
clear genome_file nAnkle nHip maxAnkle maxHip Mamp mamp N...
    Mw mw Keys Range MutDelta0 MutDelta1

% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01;
MML.tEnd = 30; % 15

%% for symmetric 2neuron CPG "knee checking" Training
% check MoE results, 'tau','b' and 'c' fixed 
% 'a' is changing within a certain range

expertCount = 10;
numOfIteretions = 10;

paper_MoE_flag = true;

if false % save data for python use
    filename = 'dataMatrixCSV_2Nsymm_08_02_2017.csv';
    csvwrite(filename,[sampl_test;freq_fromCode;freq_Matsuoka]);
end
  
switch paper_MoE_flag
    case 1 %train the MoE with the paper's method
        [ExpertsWeights, gateWeights,trainingIds,testingIds] = paper_MoE_train(sampl, targ, expertCount, numOfIteretions, 0.05, 0.9);
        [~, freq_fromNN,gateOut] = paper_MoE_test(sampl(:,testingIds), targ(:,testingIds), ExpertsWeights, gateWeights,0);
        freq_fromCode = targ(:,testingIds);
    case 0  %train the MoE with our method
        maxEphocs = 5;
        ExpertHidLayer = 1;  ExpertHidNueron = 2;
        GateHidLayer = 1;   GateHidNueron = 2;
        competetiveFlag = 3;

        [ expertsNN,gateNet,expert_i_GroupSize,gateNN_perf_vec,...
            Experts_perf_mat,Moe_perf_over_iter,emptyGroupIndecator,trainingIds,testingIds]= ...
            my_MoE_train(sampl,targ,expertCount,numOfIteretions,maxEphocs,ExpertHidLayer,ExpertHidNueron,...
                         GateHidLayer,GateHidNueron,competetiveFlag);
        [freq_fromNN,gateOut,freq_fromCode,belongToExpert,cluster_i_train_ind] =...
            my_MoE_testNet(sampl(:,testingIds),targ(:,testingIds),expertsNN,gateNet,competetiveFlag);

%             my_MoE_plotPerf(freq_fromNN(n,:),targ1,gateOut,cluster_i_train_ind,Moe_perf_over_iter,...
%                 gateNN_perf_vec,expert_i_GroupSize,Experts_perf_mat,emptyGroupIndecator,...
%                 'both',competetiveFlag);
end

%% ploting "knee" checking 2D
% figure of matsuoka estimation performance
figure;
plot(sampl_test(3,:),freq_fromCode,'g-o'); hold on
plot(sampl_test(3,:),freq_Matsuoka,'b-o');
xlabel('a');    ylabel('freq [Hz]'); grid on;
legend('real freq','our MoE');
hold off
disp('The matsuokas estimation perf:');
[~,~] = NN_perf_calc(freq_fromCode,freq_Matsuoka,1,0,'test');

% plotting a 2 experts output (freq over 'a') with color coding
figure;
plot(sampl_test(3,:),freq_fromCode,'g-o'); hold on
plot(sampl_test(3,:),freq_Matsuoka,'b-o');
g_max = vec2ind(gateOut);
colors = [1,0,0; 0,0,1]; % colors = rand(expertCount,3); % if more than 2 experts
for j=1:expertCount
    for i=1:size(sampl_test(3,:),2)
        if g_max(1,i) == j
            plot(sampl_test(3,i),freq_fromNN(n,i),'k-o','MarkerFaceColor',colors(j,:));
        end
    end
end
xlabel('a');    ylabel('freq [Hz]'); grid on;
legendNames = {'real freq','Matsuoka est','our MoE'}; %'Matsuoka estimation'
legend(legendNames);
hold off

% % plotting with many different expert num:
% figure;
% plot(vertcat(results_2N_sim(:).a),freq_fromCode,'g-o'); hold on
% plot(vertcat(results_2N_sim(:).a),freq_Matsuoka,'b-o');
% plot(vertcat(results_2N_sim(:).a),freq_fromNN(n,:));
% xlabel('a');    ylabel('freq [Hz]'); grid on;
% legendNames = {'real freq','Matsuoka estimation',cell(1,length(expertCount))};
% for j=3:(length(expertCount)+2)
%     legendNames{1,j} = ['#',num2str(expertCount(1,j-2)),' expert'];
% end
% legend(legendNames);
% title('frequency over "a" with different #experts');

disp('The MoE perf:');
[~,~] = NN_perf_calc(freq_fromCode,(freq_fromNN(n,:)),1,0,'test');

%% plotiing "knee" checking 3D
tau_test = sampl(1,testingIds);
b_test = sampl(2,1); %should be = 2.5
a_test = sampl(3,testingIds);
c_test = sampl(4,1); %should be = 5

freq_Matsuoka = zeros(1,length(testingIds));
for i=1:length(testingIds)
    a = a_test(1,i);
    tau = tau_test(1,i);
    b = b_test ;
    T = 5*tau_test(1,i);
    freq_Matsuoka(1,i) = (1/(T*2*pi))*sqrt(((tau+T)*b-(tau*a))/(tau*a)); % in [Hz]
end

figure;
scatter3(a_test,tau_test,freq_fromCode,'o'); hold on
scatter3(a_test,tau_test,freq_Matsuoka,'x');
scatter3(a_test,tau_test,freq_fromNN,'d');
xlabel('a');    ylabel('tau');    zlabel('freq [Hz]'); grid on;

legend('freq from code','Matsuoka est','MoE estimation');
title('frequency over "a" and "tau"');

disp('The MoE perf:');
[~,~] = NN_perf_calc(freq_fromCode,freq_fromNN,1,0,'test');

disp('The matsuokas estimation perf:');
[~,~] = NN_perf_calc(freq_fromCode,freq_Matsuoka,1,1,'test');
%% Train NN with 1 layer (for comparison)
net = feedforwardnet(2);
[net, tr_net] = train(net, sampl, targ);
freq_fromNN = net(sampl(:,tr_net.testInd));
disp('The normal NN perf:');

% figure;
% plot(sampl_test(3,:),freq_fromCode,'g-o'); hold on
% plot(sampl_test(3,:),net(sampl_test),'r-o');
% plot(sampl_test(3,:),freq_Matsuoka,'b-o');
% xlabel('a');    ylabel('freq [Hz]'); grid on;
% legend('real freq','NN output');
% hold off


tau_test = sampl(1,tr_net.testInd);
a_test = sampl(3,tr_net.testInd);
targ_test = targ(1,tr_net.testInd);

freq_Matsuoka_test = zeros(1,length(tr_net.testInd));
for i=1:length(tr_net.testInd)
    a = a_test(1,i);
    tau = tau_test(1,i);
    b = 2.5 ;
    T = 5*tau_test(1,i);
    freq_Matsuoka_test(1,i) = (1/(T*2*pi))*sqrt(((tau+T)*b-(tau*a))/(tau*a)); % in [Hz]
end
[~,~] = NN_perf_calc(targ_test,freq_Matsuoka_test,1,0,'test');
[~,~] = NN_perf_calc(targ_test,freq_fromNN,1,0,'test');

figure;
scatter3(a_test,tau_test,freq_fromNN,'bo'); hold on
scatter3(a_test,tau_test,targ_test,'rx');
scatter3(a_test,tau_test,freq_Matsuoka_test,'gd');
legend('NN estimation','freq from code','Matsuoka est');
xlabel('a');    ylabel('tau');    zlabel('freq [Hz]'); grid on;
% axis([1 3.5 0.24 0.26 0 1.6]);

%% Run small range
close all; clc; clear all;
% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01;
MML.tEnd = 30; % 15

N = 1000;
% CPG parameters:
tau = ones(1,N)*0.25;%(0.5-0.1).*rand(1,N) + 0.1;
b = 2.5;
c = 5;
a = (3.4-1.4).*rand(1,N) + 1.4;

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
% %     % 2neuron symmetric specific range%%
    Keys = {'\tau_r', 'beta',     'amp_2n',        '2neuron_symm_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
                  1 ,      1,          2*N,                             1,        1 ,       2*N ,            0 };
    Range = {  0.02 ,      1,         mamp,                             1,   -0.001 ,  -0.2*Mamp; % Min
               0.6  ,    8.0,         Mamp,                             6,   0.001 ,   0.2*Mamp}; % Max

    MutDelta0 = 0.04;   MutDelta1 = 0.02;
    
    save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
        'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
        'MutDelta0', 'MutDelta1', 'Keys', 'Range');
end   % define the genome file structure for matsuokaML.runSim

disp('start with the sim:');
parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
    disp(['at sim #',num2str(i)]);
    seq = [tau(1,i),b,c,0,a(1,i),0,0,0];
    [out, ~, ~] = MML.runSim(seq);
        % Prepare output:
    % Parameters
    results(i).seq = seq;
    results(i).tau = tau(1,i);
    results(i).b = b;
    results(i).c = c;
    results(i).a = a(1,i);
    results(i).x0 = out.x0;

    % Results
    results(i).periods = out.periods;
    results(i).pos_work = out.pos_work;
    results(i).neg_work = out.neg_work;
    results(i).perError1 = out.perError1;
    results(i).perOK1 = out.perOK1;
    results(i).perError2 = out.perError2;
    results(i).perOK2 = out.perOK2;
    results(i).neuronActive = out.neuronActive;
    results(i).neuronOsc = out.neuronOsc;
end 
disp('sim end...');

save('MatsRandomRes_2Neurons_change_only_a.mat','results');
