

clear all; close all; clc

% set default options
set(0,'defaultlinelinewidth',2);

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;

InFiles_names = {'VGAM_2N_symm_09_27_15_23_GA_only.mat',...
    'VGAM_2N_symm_09_28_09_05_GA_only.mat',...
    'VGAM_2N_symm_09_27_15_28_NN_classi_only.mat',...
    'VGAM_2N_symm_09_27_15_26_rescale_only.mat'};

Legends = {'GA only 1','GA only 1','NN classi','rescale'};

% the order of the parametrs in CPG Sequence:
seqOrder = {'tau' ,'b', 'c', 'NR', 'a',...
    'k_tau','k_{c}'};
% "NR" - not relevnt param 

GA_graphs = plotMOOGA4Paper(MML,InFiles_names,Legends,seqOrder);

%% collect Seqs and Fits from all GAs:
wantedGen = 1:20;

GAnum = numel(GA_graphs.Legends);   %how many GA results files we have
popSize = GA_graphs.data{1,1}.GA.Population;
seqSize = size(GA_graphs.data{1,1}.GA.Seqs,2);
genNum = length(wantedGen);

Seqs = [];
Fits = [];
for i=1:GAnum
    for j=1:genNum
        g = wantedGen(j);
        Seqs = [Seqs; GA_graphs.data{1,i}.GA.Seqs(:,:,g)];
        Fits = [Fits; GA_graphs.data{1,i}.GA.Fit(:,:,g)];
    end
end

%% get NN in/out


targ = ( Fits(:,[3]) )';  % RangeVelFit
sampl = ( Seqs(:,[1,2,3,5,6,7]) )';


% targ = ( Fits(:,[1]) )';    % VelFit
% sampl = ( Seqs(:,[1,2,3,5]) )';

%% Neural Network #1a: (CPU training)

% % create error vector
errWeigths = 1*ones(1,length(targ));

% % % OPTION: give different weights to "good" genes:
% errWeigths(1,targ>0.01) = 1*ones(1,sum(targ>0.01));

% % % OPTION: generates noise on the targets:
targ = targ + 0.01*rand(1,length(targ));


architecture = [6,20];

net = fitnet(architecture);
% net = cascadeforwardnet(architecture);

net.trainFcn = 'trainbr';
% net.divideParam.trainRatio = 0.5;
% net.divideParam.valRatio = 0.35;
% net.divideParam.testRatio = 0.15;

net.trainParam.showWindow = 1; 
net.trainParam.showCommandLine = 1;
net.trainParam.epochs = 1000;
% net.performParam.normalization = 'percent'; %It must be 'none', 'standard' or 'percent'

[net, tr] = train(net, sampl, targ,[],[],errWeigths);
[r,m,b] = regression(targ,net(sampl));
disp(['r = ',num2str(r),'    m = ',num2str(m),'    b = ',num2str(b)]);
clear r m b

% figure;
% histogram(targ,20,'Normalization','pdf'); hold on;
% histogram(net(sampl),20,'Normalization','pdf');
% legend('targets','NN outputs');

%%
rang3VelFit = ( Fits(:,[3]) )';  % RangeVelFit

% % others
% ParamNames = {'kc','c'};
% Param = ( Seqs(:,[7,5]) )';

% others
ParamNames = {'kt','tau'};
Param = ( Seqs(:,[6,1]) )';

% % % K's
% ParamNames = {'k_{tau}','k_c'};
% Param = ( Seqs(:,[6,7]) )';

Param_good = Param(:,rang3VelFit > 0.01);
Param_bad = Param(:,rang3VelFit < 0.01);

figure;
subplot(1,2,1);
histogram2(Param_good(1,:),Param_good(2,:),...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('high fits');
xlabel(ParamNames{1,1});
ylabel(ParamNames{1,2});
subplot(1,2,2);
histogram2(Param_bad(1,:),Param_bad(2,:),...
    'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','pdf');
title('low fits');
xlabel(ParamNames{1,1});
ylabel(ParamNames{1,2});