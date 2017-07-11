
% V example:
% From: https://goker.wordpress.com/2011/07/01/mixture-of-experts/

% Generate data for mixture of experts regression test
% data can be learned with a 2 logistic discrimination expert MoE 

clear all; close all; clc
%% regression dataset
x2 = rand(5000,1);
r2 = zeros(5000,1);
r2(x2<0.5,:) = x2(x2<0.5) .* 2 + ( randn(size(x2(x2<0.5))) .* 0.1 );
r2(x2>=0.5,:) = (2 - (x2(x2>=0.5) .* 2)) + ( randn(size(x2(x2>=0.5))) * 0.1 );
tx2 = x2(1:3000,:);
vx2 = x2(3001:500,:);
tr2 = r2(1:3000,:);
vr2 = r2(3001:500,:);

% plot(x2, r2, 'ob');

sampl=x2';
targ=r2';

clear x2 r2 tx2 vx2  tr2 vr2

%% using the paper script
expertCount = 2;
numOfIteretions = 30;

[ExpertsWeights, gateWeights,errs] = paper_MoE_train(sampl, targ, expertCount, numOfIteretions,0.1, 0.99995);

[~, freq_fromNN,g] = paper_MoE_test(sampl,targ, ExpertsWeights, gateWeights,1);

paper_MoE_plotPerf(sampl,targ,freq_fromNN,g,errs);



%% My code:
expertCount = 2;      % how many "experts" (fitting NN)
numOfInputs = 1; %how many inputs to each expert
maxEphocs = 5;      % max number of ephocs for each NN training
numOfIteretions = 10;  % number of loop interations
ExpertHidLayer = 1; % num of hidden layer in each expert
ExpertHidNueron = 2; % num of neurons in each hidden layer
GateHidLayer = 1; % num of hidden layer in gateNN
GateHidNueron = 2; % num of neurons in each hidden layer
competetiveFlag = 1; % if '1'- "winner takes all"
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


if competetiveFlag==1 || competetiveFlag==2
    
    [~,~] = NN_perf_calc(targ,netOut,1,0);
    samplNew = [];
    for j=1:expertCount % rearenging the samples to fit the targets
        sampl_temp = sampl(:,cluster_i_train_ind{1,j});
        samplNew = [samplNew,sampl_temp];
    end
    sampl = samplNew;
    
    
    colors = [1,0,0;
              0,0,1];
    markers = ['o','x'];
    legendNames = cell(1,expertCount);
    for j=1:expertCount
        legendNames{1,j} = ['#',num2str(j),' expert'];
    end
    figure; hold on
    for j=1:expertCount
        out_temp = netOut(:,cluster_i_train_ind{1,j});
        sampl_temp = sampl(:,cluster_i_train_ind{1,j});

        h = plot(sampl_temp,out_temp,'Color',colors(j,:),'LineStyle','none');
        h.Marker = markers(1,j);
    end
    hold off;
    legend(legendNames);
    title('which sample belongs to which expert');
    
end

figure; 
h2 = plot(sampl,targ,'LineStyle','none'); hold on;
h2.Marker = 'o';
h3 = plot(sampl,netOut,'Color','r','LineStyle','none');
h3.Marker = 'x';
legend('targets','output');


