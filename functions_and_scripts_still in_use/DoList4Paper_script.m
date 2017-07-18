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

% figure;
% h=histogram(periods_ratios,100); grid minor;
% h.BinLimits = [0,2.5];
% h.BinWidth = 0.1;
% h.Normalization = 'pdf';
% xlabel('P_{hip} / P_{ankle} : the ratio between the two CPG outputs');
% ylabel('probability density');
% title('histogram of the ratio between the two CPG outputs');
% set(gca,'FontSize',10);
% savefig('figure_TBD_Histogram_of_ratio_between_periods_hipAnkle')


ids = osc_ids & diff_ids;

periodsMean = mean(periods(ids),1);
%% Figure 5:
clc

results_osc = results(ids);
results_n_osc = results(~ids);

seq_osc = (vertcat(results_osc(:).seq))';
seq_osc = seq_osc(1:18,:);

seq_n_osc = (vertcat(results_n_osc(:).seq))';
seq_n_osc = seq_n_osc(1:18,:);

for i = 1:length(seqOrder)
    p_name = seqOrder{1,i};
    
    p_vec = seq_n_osc(strcmp(p_name,seqOrder),:);
    p_osc_vec = seq_osc(strcmp(p_name,seqOrder),:);
    
    hist_compare(p_vec,p_osc_vec,p_name,20,{'n-osc CPGs','osc CPGs'},'plot');

end

clear p_name p_vec p_osc_vec pvalue rejection

%% "Kullback-Leibler" divergence:
clc


p_name = 'tau'; % seqOrder{1,i};

p_vec = sort(seq_n_osc(strcmp(p_name,seqOrder),:));
p_osc_vec = sort(seq_osc(strcmp(p_name,seqOrder),:));
    
temp = ();


clear p_name p_vec p_osc_vec pvalue rejection

%% appendix to Figure 5
% join all 'C_i' to one vector and all 'W_ij' to another vector and 
%   compare the distribution (instead of using many Hstograms).

% seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
%     'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
%     'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};

% join C_i:
c_i = [];
c_i_osc = [];

for i = 3:6
    p_name = seqOrder{1,i};
    
    p_vec = seq_n_osc(strcmp(p_name,seqOrder),:);
    p_osc_vec = seq_osc(strcmp(p_name,seqOrder),:);

    c_i = [c_i,p_vec];
    c_i_osc = [c_i_osc,p_osc_vec];
    
    clear p_name p_vec _osc_vec
end
hist_compare(c_i,c_i_osc,'all c_{i}',20,{'n-osc CPGs','osc CPGs'},'plot');
clear c_i c_i_osc

% join W_ij:
W_ij = [];
W_ij_osc = [];
for i = 7:18
    p_name = seqOrder{1,i};
    
    p_vec = seq_n_osc(strcmp(p_name,seqOrder),:);
    p_osc_vec = seq_osc(strcmp(p_name,seqOrder),:);

    W_ij = [W_ij,p_vec];
    W_ij_osc = [W_ij_osc,p_osc_vec];
    
    clear p_name p_vec _osc_vec
end
hist_compare(W_ij,W_ij_osc,'all original W_{ij}',20,{'n-osc CPGs','osc CPGs'},'plot');
kullback_leibler_divergence(W_ij',W_ij_osc')
clear W_ij W_ij_osc

%% check normalized Matsuoka coupling weights:
Wnames = {'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};
CNames = {'c_1','c_2','c_3','c_4'};

% for oscillating CPG's:
[Worig_osc,c_osc] = prepare_NN_inOut(seq_osc,periodsMean,Wnames,CNames,seqOrder);
W_norm_osc = Worig_to_W_normalize(Worig_osc,c_osc );

% for non-oscillating CPG's:
[Worig_n_osc,c_n_osc] = prepare_NN_inOut(seq_n_osc,periodsMean,Wnames,CNames,seqOrder);
W_norm_n_osc = Worig_to_W_normalize(Worig_n_osc,c_n_osc );

W_ij = [];
W_ij_osc = [];
for i = 1:12
    p_vec = W_norm_osc(i,:);
    p_osc_vec = W_norm_n_osc(i,:);

    W_ij = [W_ij,p_vec];
    W_ij_osc = [W_ij_osc,p_osc_vec];
    
    clear p_name p_vec _osc_vec
end
hist_compare(W_ij,W_ij_osc,'all normalized W {ij}',20,{'n-osc CPGs','osc CPGs'},'plot');

clear W_ij W_ij_osc
clear Wnames CNames Worig_osc c_osc Worig_n_osc c_n_osc W_norm_osc W_norm_n_osc
%% Figure 6:

inputsNames = {'periods','tau'...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};
outputsNames = {'b'};

NumOfRepeats = 5;

% sampl = zeros(length(inputsNames),size(seq_osc,2));
% targ = zeros(length(outputsNames),size(seq_osc,2));

[sampl,targ] = prepare_NN_inOut(seq_osc,periodsMean,inputsNames,outputsNames,...
    seqOrder);

% perf over samples num:
outputFileName = 'NN_Perf_over_sampl_num_test.mat';
HiddenN = 10;
numOfSampl = 1000*[20,30,50,75,100,150,200];
out_sampl_num = NN_Perf_over_samples_num(sampl,targ,NumOfRepeats,HiddenN,...
    numOfSampl,'train',outputFileName);
NN_Perf_over_samples_num(sampl,targ,NumOfRepeats,HiddenN,...
    numOfSampl,'plot',outputFileName);
NN_Perf_over_samples_num(sampl,targ,NumOfRepeats,HiddenN,...
    numOfSampl,'text',outputFileName);

% perf over hid neuron num:
outputFileName = 'NN_Perf_over_HNnum_test_1.mat';
HiddenN = [2,4,6,8,10,20,30,40,50,60,70,80];
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

clc 

inputsNames_4GA = {'period_desired','tau'...
    'w_{12}','w_{13}','w_{14}',...
    'w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}',...
    'w_{41}','w_{42}','w_{43}'};
outputsNames_4GA = {'b'};

% find NOT oscilatory CPGs:
not_osc_ids = find(~osc_ids); 
if length(not_osc_ids)<500
    results_n_osc = load('MatsRandomRes_4Neurons_with_LSQ_1_2.mat','results');
    results_n_osc = results_n_osc.results;
    periods = horzcat(results_n_osc(:).periods);
    non_osc_ids = isnan(periods);
    non_osc_ids = non_osc_ids(1,:) & non_osc_ids(2,:);
    cpg_non_osc = randsample(find(non_osc_ids),500);
    results_old = results_n_osc(cpg_non_osc);
    clear results_n_osc
else
    cpg_non_osc = randsample(not_osc_ids,500); % use 500 n-osc CPGs
    results_old = results(cpg_non_osc);
end


seq_n_osc = (vertcat(results_old (:).seq))';
seq_n_osc = seq_n_osc(1:18,:);


[sampl_4GA,targ_4GA] = prepare_NN_inOut(seq_n_osc,cpg_non_osc,...
    inputsNames_4GA,outputsNames_4GA,seqOrder);

caseNum = 7;
% HiddenN = [10,20,30,40,50];
HiddenN = [2,5,10,15];
samplNum = size(sampl,2);
trainRatio = 0.7;
valRatio = 0.15;
testRatio = 1-trainRatio-valRatio;

numRepeat = 10;

accuracy = zeros(numRepeat,length(HiddenN));
percent_osc_new = zeros(numRepeat,length(HiddenN));
conv_in_range = zeros(numRepeat,length(HiddenN));
 
for i=1:length(HiddenN)
    disp(['NN with ',num2str(HiddenN(1,i)),' hidden neurons:']);
    for j=1:numRepeat
        disp(['    iter ',num2str(j),' out of ',num2str(numRepeat)]);
        % Shuffle the samples (matrix rows):
        [trainInd,valInd,testInd] = dividerand(samplNum,trainRatio,valRatio,testRatio);


        % prepare the Neural Network:
        net = feedforwardnet(HiddenN(1,i));

        % Neural Network training parameters:
        net.trainParam.showWindow = false; % dont show training window
        net.divideFcn = 'divideind';
        net.divideParam.trainInd = trainInd;
        net.divideParam.valInd   = valInd;
        net.divideParam.testInd  = testInd;

        % NN training
        if 1
        net = train(net, sampl, targ);
        else % use untrained net (with random weights)
        net = configure(net,sampl,targ);
        end
        
        [percent_osc_new(j,i),conv_in_range(j,i),accuracy(j,i)] = ...
             NN_GA_perf(net,sampl_4GA,results_old,seqOrder,caseNum);
    end
end

percent_osc_new_mean = mean(percent_osc_new,1);
conv_in_range_mean = mean(conv_in_range,1);
accuracy_mean = mean(accuracy,1);

percent_osc_new_std = std(percent_osc_new,[],1);
conv_in_range_std = std(conv_in_range,[],1);
accuracy_std = std(accuracy,[],1);

means = [percent_osc_new_mean;...
    conv_in_range_mean;...
    accuracy_mean];
stdevs = [percent_osc_new_std;...
    conv_in_range_std;...
    accuracy_std];
% Names = {' ','10neurons',' ','20neurons ',' ','30neurons ',' ','40neurons ',' ','50neurons'};
Names = {' ','2neurons',' ','5neurons ',' ','10neurons ',' ','15neurons '};
label_Y = '';
graph_title = ['NN perf over different hidden neurons num'];
graph_legend = {'converged','conv in range','accuracy'};
plot_bars_with_errors(means,stdevs,...
    Names,label_Y,graph_title,graph_legend);

disp('|-------------------------------------------------------------|');
disp('| hidNnum | conv mean | conv std | inRange mean | inRange std |');
disp('|-------------------------------------------------------------|');
res = cell(1,length(HiddenN));
for i=1:length(HiddenN)
    res{1,i} = sprintf('|case %d   | %0.3f     | %0.3f    | %0.3f        | %0.3f       | \n',HiddenN(1,i),...
        means(1,i),stdevs(1,i),means(2,i),stdevs(2,i));
    disp(res{1,i});
end
disp('|-------------------------------------------------------------|');