clear all; close all; clc


%% Load Data
load('MatsRandomRes_2Neurons_general_16_01_2017_A.mat','results')
results1 = results; periods1 = horzcat(results(:).periods);
clear results
load('MatsRandomRes_2Neurons_general_16_01_2017_B.mat','results')
results2 = results; periods2 = horzcat(results(:).periods);
clear results
load('MatsRandomRes_2Neurons_general_16_01_2017_C.mat','results')
results3 = results; periods3 = horzcat(results(:).periods);
clear results
load('MatsRandomRes_2Neurons_general_16_01_2017_D_narrowRange.mat','results')
results4 = results; periods4 = horzcat(results(:).periods);
clear results
% concetrate the data:
results = horzcat(results1,results2,results3);%,results4);
periods = horzcat(periods1,periods2,periods3);%,periods4);

clear results1 results2 results3% results4
clear periods1 periods2 periods3% periods4

%% identify which neurons have periods
ids_period = ~isnan(periods); % only ones with period
ids_error = (max(horzcat(results(:).perError2)',[],2) < 0.001)'; % only ones with low enought error
ids = find(ids_period & ids_error);

%% prepare NN inputs and outputs:
% parametersCells = {'tau','b','a','s'};
% % parametersCells = {'tau','b'};
% targetCells = {'freq'};
% [ sampl,targ ] = prepareData_2Neurons( results,periods,ids,parametersCells,targetCells);
% %
% figure;
% histogram(max(horzcat(results(ids).perError2)',[],2),1000);
% title('Hist of period error'); xlabel('period error');

parametersCells = {'tau','b','w12','w21'};
targetCells = {'freq'};
sampl = vertcat(results(ids).seq)'; % extracting the parameters from the structure to a matrix
sampl = sampl([1,2,5,6],:); % get rid of the the unimportant 'gen' (the k's for the control).
targ= 1./periods(ids);

clear ids_period ids_error nSims
% clear results periods

%% save to xlse file for eureqa:
filename = 'dataMatrixXLSX_2neurons_updated.xlsx';
dataMatrix = [sampl;targ];
dataMatrix = dataMatrix(:,1:10000);
xlswrite(filename,dataMatrix');
%% normalized inputs
normParams = zeros(size(sampl, 1), 2);
notNormSampl = sampl;
for i = 1:size(sampl, 1)
    feat = sampl(i, :);
    normParams(i, :) = [mean(feat), std(feat)];
    sampl(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
    
end

%% Performance over Hidden neuron Number
NumOfRepeats = 5;
HiddenN = [2,3,4,5,6,7,8,9,10,11,12];
[ NN_Mean_over_HN_num,NN_stdev_over_HN_num ] = NN_Perf_over_HNnum( NumOfRepeats,HiddenN,sampl,targ,1 );

%% Performance over samples quantity
NumOfRepeats = 7;
HiddenN = 10;
dataPointsNum = [1000,5000,7000,10000,15000,20000,30000,50000,70000,...
    100000];
[ NN_Mean_over_sampl_num,NN_stdev_over_sampl_num ] = NN_Perf_over_sampl_num( NumOfRepeats,HiddenN,dataPointsNum,sampl,targ,1 );


%% visualiesd the results in the W matrix
NN_weights_matrix_plot( net,parametersCells )

%% Train NN with 1 layer
HiddenN = 5;
net = feedforwardnet(HiddenN);
% net = fitnet(HiddenN);
[net, tr] = train(net, sampl, targ);

%% Train NN with 2 layers
HiddenN = 10;
net = feedforwardnet([HiddenN,HiddenN]);
[net, tr] = train(net, sampl, targ);
%% Train NN with 3 layers
HiddenN = 3;
net = feedforwardnet([HiddenN,HiddenN,HiddenN]);
[net, tr] = train(net, sampl, targ);

%% Calculating the R^2:

inputs = sampl;
outputs = targ;
netOut = net(inputs);
err = outputs-netOut;
errVar = var(err,0,2);
inputVar = var(outputs,0,2);

R_squar = 1-(errVar/inputVar);

clear inputs outputs netOut err errVar inputVar

%%
%% run many genes with similar tau,T,b,c and different a and check the period
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
    % %     % Final genome with tau_r + beta (constant tau_u/tau_v ratio) %%
        Keys = {'\tau_r', 'beta', 'amp',        '2neuron_general_weights', 'ks_\tau',     'ks_c', 'IC_matsuoka';
                      1 ,      1,   2*N,                                2,        1 ,       2*N ,            0 };
        Range = {  0.02 ,    0.2,  mamp,                            [1,1],   -0.001 ,  -0.2*Mamp; % Min
                   0.6  ,   10.0,  Mamp,                          [10,10],   0.001 ,   0.2*Mamp}; % Max

    MutDelta0 = 0.04;   MutDelta1 = 0.02;

    save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
        'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
        'MutDelta0', 'MutDelta1', 'Keys', 'Range');
end % define the genome structure
   
close all; clc; %clear all;
MML = MatsuokaML(); % calling Matsuoka class
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01;
MML.tEnd = 30; % 15

N = 25;
freq_fromNN = zeros(1,N*N);
freq_fromCode = zeros(1,N*N);
tau = 0.25;
b = 2.5;
c = 5;
w12 = linspace(1.6,3.4,N);
w21 = linspace(1.6,3.4,N);
clear results peiords
tic
for i=1:N
    for j=1:N
        seq = [tau,b,c,c,w12(1,i),w21(1,j),0,0,0];
        [out, ~, ~] = MML.runSim(seq);
            % Prepare output:
        % Parameters
        k = j+N*(i-1);
        results(k).seq = seq;
        results(k).tau = tau;
        results(k).b = b;
        results(k).c = c;
        results(k).w12 = w12(1,i);
        results(k).w21 = w21(1,j);
        results(k).x0 = out.x0;

        % Results
        results(k).periods = out.periods;
        results(k).pos_work = out.pos_work;
        results(k).neg_work = out.neg_work;
        results(k).perError1 = out.perError1;
        results(k).perOK1 = out.perOK1;
        results(k).perError2 = out.perError2;
        results(k).perOK2 = out.perOK2;
        results(k).neuronActive = out.neuronActive;
        results(k).neuronOsc = out.neuronOsc;
        
        freq_fromNN(1,k) = net([tau;b;w12(1,i);w21(1,j)]);
        freq_fromCode(1,k)= 1/out.periods;
        
        clear out
    end
end
toc
%
w12_vec = vertcat(results(:).w12);
w21_vec = vertcat(results(:).w21);
figure;
scatter3(w12_vec,w21_vec,freq_fromNN,'b','d','filled'); hold on;
scatter3(w12_vec,w21_vec,freq_fromCode,'k','o','filled');
xlabel('w_{12}');    ylabel('w_{21}');    zlabel('freq [Hz]'); grid on;
legend('blue: freq from NN','black: freq from Code');
title('frequency over w_{12} and w_{21}');

vv = randsample(w12,4); vv=vv';
figure;
for i=1:4
    indx = find(w12_vec==vv(i));
    subplot(2,2,i);
    scatter(vertcat(results(indx).w21),freq_fromNN(indx),'b'); hold on
    scatter(vertcat(results(indx).w21),freq_fromCode(indx),'r');
    xlabel('w_{21}');   ylabel('freq[Hz');  title(['w_{12} = ',num2str(vv(i))]);
    
end
legend('blue: freq from NN','red: freq from Code');

figure;
indx2 = find(w12_vec==w21_vec);
wii_vec = w12_vec(indx2);
scatter(wii_vec,freq_fromNN(indx2),'b'); hold on
scatter(wii_vec,freq_fromCode(indx2),'r');
xlabel('w_{12} = w_{21}');   ylabel('freq[Hz');  title('w_{12} = w_{21}');
legend('blue: freq from NN','red: freq from Code');
