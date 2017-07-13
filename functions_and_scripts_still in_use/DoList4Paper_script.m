% 
% compare the probability density of each CPG parameters (both osc and
% non-osc).
% 
% 
% 
clear all; close all; clc

load('MatsRandomRes_4Neurons_with_LSQ_1_2.mat','results');
% load('MatsRandomRes_all_from_1-2_2017.mat','results');

seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};

periods = horzcat(results(:).periods);

%% get only good CPG's:

% Filter CPG's where not both signals oscillating:
osc_ids = ~isnan(periods);
osc_ids = osc_ids(1,:) & osc_ids(2,:);

% Filter CPG's where the is a big difference between hip and ankle:
periods_ratios = (periods(1,:)./periods(2,:));
diff_ids = (periods_ratios >  0.85) & (periods_ratios <  1.15); 

figure;
h=histogram(periods_ratios,100); grid minor;
h.BinLimits = [0,2.5];
h.BinWidth = 0.1;
h.Normalization = 'pdf';
xlabel('P_{hip} / P_{ankle} : the ratio between the two CPG outputs');
ylabel('probability density');
title('histogram of the ratio between the two CPG outputs');
set(gca,'FontSize',10);

figure;
[Est_pdf, Edges_pdf] = ...
    histcounts(periods_ratios,1000, 'Normalization','pdf');
bin_center = (Edges_pdf(1:end-1)+ Edges_pdf(2:end))/2;
bar(bin_center,Est_pdf);
clear Est_pdf Edges_pdf bin_center


ids = osc_ids & diff_ids;

periodsMean = mean(periods(ids),1);
%% Figure 5:
results_osc = results(ids);
seq_osc = (vertcat(results_osc(:).seq))';
seq_osc = seq_osc(1:18,:);

seq = (vertcat(results(:).seq))';
seq = seq(1:18,:);

for i = 1:length(seqOrder)
    p_name = seqOrder{1,i};
    
    p_vec = seq(strcmp(p_name,seqOrder),:);
    p_osc_vec = seq_osc(strcmp(p_name,seqOrder),:);
    
    [pvalue,rejection] = ranksum(p_vec,p_osc_vec);
    str1 = sprintf('p_value of %s is %0.3f and the hypotesis got %d',...
        p_name,pvalue,rejection);
    disp([str1,' real pval is ',num2str(pvalue)]);
    
%     figure;
%     histogram(p_vec,100,'Normalization','pdf'); grid minor; hold on;
%     histogram(p_osc_vec,100,'Normalization','pdf');
%     xlabel(p_name);
%     ylabel('probability density');
%     title(['histogram of ',p_name]);
%     legend('all CPGs','Osc CPGs');
%     set(gca,'FontSize',13);
%     hold off;
%       
end

clear p_name p_vec p_osc_vec pvalue rejection

%% Figure 6:

inputsNames = {'period_desired','tau'...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};
outputsNames = {'b'};

NumOfRepeats = 5;

sampl = zeros(length(inputsNames),size(seq_osc,2));
targ = zeros(length(outputsNames),size(seq_osc,2));

% Prepare inputs:
for i = 1:length(inputsNames)
    p_name = inputsNames{1,i};
    switch p_name
        case 'period_desired'
            per_des_Min = 0.68;
            per_des_Max = 0.78;
            sampl(i,:) = per_des_Min + ...
                ((per_des_Max-per_des_Min) * rand(1,sum(ids)));
            clear per_des_Min per_des_Max
        otherwise
            sampl(i,:) = seq_osc(strcmp(p_name,seqOrder),:);
    end   
end
% Prepare inputs and outputs:
for i = 1:length(outputsNames)
    p_name = outputsNames{1,i};
    
    targ(i,:) = seq_osc(strcmp(p_name,seqOrder),:);
    clear p_name
end

% perf over samples num:
outputFileName = 'NN_Perf_over_sampl_num_test.mat';
HiddenN = 5;
numOfSampl = 1000*[10,20,30];
out_sampl_num = NN_Perf_over_samples_num(sampl,targ,NumOfRepeats,HiddenN,...
    numOfSampl,'train',outputFileName);
NN_Perf_over_samples_num(sampl,targ,NumOfRepeats,HiddenN,...
    numOfSampl,'plot',outputFileName);
NN_Perf_over_samples_num(sampl,targ,NumOfRepeats,HiddenN,...
    numOfSampl,'text',outputFileName);

% perf over hid neuron num:
outputFileName = 'NN_Perf_over_HNnum_test.mat';
HiddenN = [2,3,4];
keepRatioConstant = true;
out_HidN_num = NN_Perf_over_Hid_Neuron_Num(sampl,targ,NumOfRepeats,...
    HiddenN,'train',keepRatioConstant,outputFileName);
NN_Perf_over_Hid_Neuron_Num(sampl,targ,NumOfRepeats,...
    HiddenN,'plot',keepRatioConstant,outputFileName);
NN_Perf_over_Hid_Neuron_Num(sampl,targ,NumOfRepeats,...
    HiddenN,'text',keepRatioConstant,outputFileName);

%% Figure 8:
% non-osc CPGs that produced oscillations after pre-tuning with NNs
% (of selected configuration, currently configuration #4) 
% with different hidden neurons.
% THIS SHOULD JUSTIFY THE SELECTION OF THE NN FOR MOGA.


% find NOT oscilatory CPGs:
not_osc_ids = find(~osc_ids); 
cpg_non_osc = randsample(not_osc_ids,500); % use 500 n-osc CPGs
results_old = results(cpg_non_osc);

caseNum = 7;
HiddenN = [2,3,4];
samplNum = size(sampl,2);
trainRatio = 0.7;
valRatio = 0.15;
testRatio = 1-trainRatio-valRatio;

numRepeat = 5;

accuracy = zeros(numRepeat,length(HiddenN));
percent_osc_new = zeros(numRepeat,length(HiddenN));
conv_in_range = zeros(numRepeat,length(HiddenN));
 
for i=1:length(HiddenN)
    for j=1:numRepeat
        % Shuffle the samples (matrix rows):
        [trainInd,valInd,testInd] = dividerand(samplNum,trainRatio,valRatio,testRatio);


        % prepare the Neural Network:
        net = feedforwardnet(HiddenN);

        % Neural Network training parameters:
        net.trainParam.showWindow = false; % dont show training window
        net.divideFcn = 'divideind';
        net.divideParam.trainInd = trainInd;
        net.divideParam.valInd   = valInd;
        net.divideParam.testInd  = testInd;

        % NN training
        net = train(net, sampl, targ);

        [sampl4change,~] = prepare_NN_inOut(seq,periodsMean,inputsNames,...
            outputsNames,seqOrder);
        [percent_osc_new(j,i),conv_in_range(j,i),accuracy(j,i)] = ...
             NN_GA_perf(net,sampl4change,results_old,seqOrder,caseNum);
    end
end

