% 
% compare the probability density of each CPG parameters (both osc and
% non-osc).
% 
% 
% 
clear all; close all; clc

% load('MatsRandomRes_4Neurons_4Paper.mat','results');
% load('MatsRandomRes_4Neurons_with_LSQ_1_2.mat','results');
% load('MatsRandomRes_all_from_1-2_2017.mat','results');

% % % % Load peridos that mainly oscillating in range:
% load('MatsRandomRes_4Neurons_4Paper_Rescaled_Sims.mat','results_rescaled');
% results1 = results_rescaled;    clear results_rescaled
% load('MatsRandomRes_4Neurons_4Paper_Rescaled_Sims_2.mat','results_rescaled');
% results2 = results_rescaled';    clear results_rescaled
% results = [results1,results2];  clear results1 results2

% % Load peridos oscillating and periods which oscillates in range:
load('MatsRandomRes_4Neurons_4Paper_Rescaled_Sims.mat','results_rescaled');
results1 = results_rescaled;    clear results_rescaled
load('MatsRandomRes_4Neurons_4Paper_Rescaled_Sims_2.mat','results_rescaled');
results2 = results_rescaled';    clear results_rescaled
load('MatsRandomRes_4Neurons_4Paper.mat','results');
results3 = results;    clear results
results = [results1,results2,results3];  clear results1 results2 results3

% % % Load peridos that osc and osc in range:
% % % % for plotting histograms in range and out of rang
% load('MatsRandomRes_4Neurons_4Paper_Rescaled_Sims.mat','results_rescaled');
% results1 = results_rescaled;    clear results_rescaled
% load('MatsRandomRes_4Neurons_4Paper.mat','results');
% results2 = results;    clear results
% results = [results1,results2];  clear results1 results2

seqOrder = {'tau','b','c_1','c_2','c_3','c_4',...
    'w_{12}','w_{13}','w_{14}','w_{21}','w_{23}','w_{24}',...
    'w_{31}','w_{32}','w_{34}','w_{41}','w_{42}','w_{43}'};

periods = horzcat(results(:).periods);

% define the class for CPG simulation:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 4;
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

periodsMean = mean(periods(:,ids),1);

ids_des_period = ids & ((periods(1,:) > 0.6) & (periods(1,:) < 0.86));

ids_not_des_period = ids & ((periods(1,:) < 0.6) | (periods(1,:) > 0.86));

results_osc = results(ids);     % oscillatory CPGs
results_n_osc = results(~ids);  % non-oscillatory CPGs
results_in_range = results(ids_des_period); % oscillatory CPGs in period range
results_not_in_range = results(ids_not_des_period);

seq_osc = (vertcat(results_osc(:).seq))';
seq_osc = seq_osc(1:18,:);

seq_n_osc = (vertcat(results_n_osc(:).seq))';
seq_n_osc = seq_n_osc(1:18,:);

seq_in_range = (vertcat(results_in_range(:).seq))';
seq_in_range = seq_in_range(1:18,:);
periods_in_range = mean(periods(ids_des_period),1);

seq_not_in_range = (vertcat(results_not_in_range(:).seq))';
seq_not_in_range = seq_not_in_range(1:18,:);


%% Figure 5: (and "Kullback-Leibler" divergence)
clc

for i = 1:length(seqOrder)
    p_name = seqOrder{1,i};
    
    p_vec = seq_n_osc(strcmp(p_name,seqOrder),:);
    p_osc_vec = seq_osc(strcmp(p_name,seqOrder),:);
    
    hist_compare(p_vec,p_osc_vec,p_name,20,{'n-osc CPGs','osc CPGs'},'plot');

end

clear p_name p_vec p_osc_vec pvalue rejection

% % % % % "Kullback-Leibler" divergence:
clc

tau_n_osc = seq_n_osc(strcmp('tau',seqOrder),:);
tau_osc = seq_osc(strcmp('tau',seqOrder),:);
dist_tau = KL_div_4paper(tau_n_osc,tau_osc)
clear tau_n_osc tau_osc

b_n_osc = seq_n_osc(strcmp('b',seqOrder),:);
b_osc = seq_osc(strcmp('b',seqOrder),:);
dist_b = KL_div_4paper(b_n_osc,b_osc)
clear b_n_osc b_osc
%% appendix to Figure 5
% join all 'C_i' to one vector and all 'W_ij' to another vector and 
%   compare the distribution (instead of using many Hstograms).

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
dist_c_i = KL_div_4paper(c_i,c_i_osc)
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
dist_W_ij = KL_div_4paper(W_ij,W_ij_osc)
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

% truncted the distributions:
ids_W_ij = (W_ij < 100);    ids_W_ij_osc = (W_ij_osc < 100);

disp('the number of CPGs with abnormaly large weights after normalization:');
disp(['non oscillating: ',num2str(sum(~ids_W_ij_osc))]);
disp(['oscillating: ',num2str(sum(~ids_W_ij))]);

hist_compare(W_ij(ids_W_ij),W_ij_osc(ids_W_ij_osc),...
    'all normalized W {ij}',20,{'n-osc CPGs','osc CPGs'},'plot');
dist_W_ij_norm = KL_div_4paper(W_ij(ids_W_ij),W_ij_osc(ids_W_ij_osc))

clear W_ij W_ij_osc ids_W_ij ids_W_ij_osc
clear Wnames CNames Worig_osc c_osc Worig_n_osc c_n_osc W_norm_osc W_norm_n_osc

%% heat map heatogram for parameter 'b':
close all; clc

p_n_osc = seq_n_osc(strcmp('b',seqOrder),:);
p_osc = seq_osc(strcmp('b',seqOrder),:);
Title = ['parameter: ','b',' of n-osc CPGs'];

Edges = [];

figure;
ax1 = subplot(2,1,1);
[ax] = heatmap_histogram(ax1,p_n_osc,Edges,Title);
ax2 = subplot(2,1,2);
[ax] = heatmap_histogram(ax2,p_osc,Edges,Title);

%% Figure 4:
% this data matched to file: 'MatsRandomRes_4Neurons_4Paper.mat'

res = load('MatsRandomRes_4Neurons_4Paper_Rescaled_Sims.mat');
results_before_rescale = res.results_osc;
results_after_rescale = res.results_rescaled;
clear res


CPG_num = size(periods,2);
disp(['total number of CPGs is: ',num2str(CPG_num)]);
disp(['total number of osc CPGs is: ',num2str(size(seq_osc,2)),...
    ' which is: ',num2str(100*size(seq_osc,2)/CPG_num),'%']);
disp(['total number of CPGs with desired periods is: ',...
    num2str(sum(ids_des_period)),' which is: ',...
    num2str(100*sum(ids_des_period)/CPG_num),'%']);


% plot distribution in heatmap:
periods_before_rescale = horzcat(results_before_rescale(:).periods);
periods_after_rescale = horzcat(results_after_rescale(:).periods);

ids1 = (periods_before_rescale(1,:) > 0.3) &...
    (periods_before_rescale(1,:) < 7);
vec1 = periods_before_rescale(1,ids1);

ids2 = (periods_after_rescale(1,:) > 0.3) &...
    (periods_after_rescale(1,:) < 7);
vec2 = periods_after_rescale(1,ids2);

Title1 = sprintf('Period distribution of random parameters');
Title2 = sprintf('Period distribution of random parameters after re-scaling');

Edges = 0:0.1:7;

% figure;
% ax1 = subplot(2,1,1);
% [ax] = heatmap_histogram(ax1,vec1,Edges,Title1);
% ax2 = subplot(2,1,2);
% [ax] = heatmap_histogram(ax2,vec2,Edges,Title1);

hist_compare(vec2,vec1,'period',...
    50,{'after rescale','all osc'},'plot')

clear ax1 ax2 ax vec1 vec2 Title1 Title2 
clear ids1 ids2 periods_b4_rescale periods_rescale
%% plot distribution in heatmap:

res = results_osc(1:1000);

periods_b4_rescale = horzcat(res(:).periods);
periods_rescale = horzcat(results_rescaled(:).periods);

vec1 = periods_b4_rescale(1,:);
vec2 = periods_rescale(1,:);

Title1 = sprintf('Period distribution of random parameters');
Title2 = sprintf('Period distribution of random parameters after re-scaling');

figure;
ax1 = subplot(2,1,1);
[ax] = heatmap_histogram(ax1,vec1,100,Title1);
ax2 = subplot(2,1,2);
[ax] = heatmap_histogram(ax2,vec2,100,Title2);

clear ax1 ax2 ax vec1 vec2 periods_b4_rescale periods_rescale

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

caseNum = 9; % caseNum = 7;
% HiddenN = [10,20,30,40,50];
HiddenN = [10,20];
trainRatio = 0.7;
valRatio = 0.15;
testRatio = 1-trainRatio-valRatio;


[inputsNames_4NN,outputsNames_4NN] = ...
    check_NN_case_for_paper(caseNum,'period');

[sampl,targ] = prepare_NN_inOut(seq_in_range,periods_in_range,...
    inputsNames_4NN,outputsNames_4NN,seqOrder);

[inputsNames_4GA,outputsNames_4GA] = ...
    check_NN_case_for_paper(caseNum,'period_desired');

% number of training samples:
samplNum = size(sampl,2);

% find NOT oscilatory CPGs:
load('MatsRandomRes_4Neurons_4Paper_not_osc_Sims.mat','results_not_osc');
howMuch = length(results_not_osc);
cpg_non_osc = randsample(1:howMuch,500);
results_old = results_not_osc(cpg_non_osc);
clear results_not_osc howMuch

seq_n_osc_4Sim = (vertcat(results_old (:).seq))';
seq_n_osc_4Sim = seq_n_osc_4Sim(1:18,:);


[sampl_4GA,targ_4GA] = prepare_NN_inOut(seq_n_osc_4Sim,cpg_non_osc,...
    inputsNames_4GA,outputsNames_4GA,seqOrder);

numRepeat = 5;

accuracy1 = zeros(numRepeat,length(HiddenN));
accuracy2 = zeros(numRepeat,length(HiddenN));
accuracy3 = zeros(numRepeat,length(HiddenN));
percent_osc_new = zeros(numRepeat,length(HiddenN));
conv_in_range = zeros(numRepeat,length(HiddenN));
MSE_testErr = zeros(numRepeat,length(HiddenN));

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
            [net,tr] = train(net, sampl, targ);
            MSE_testErr(j,i) = tr.best_tperf;
        else % use untrained net (with random weights)
            net = configure(net,sampl,targ);
            MSE_testErr(j,i) = NaN;
        end
        
        [percent_osc_new(j,i),conv_in_range(j,i),...
            accuracy1(j,i),accuracy2(j,i),accuracy3(j,i)] = ...
             NN_GA_perf(MML,net,sampl_4GA,results_old,seqOrder,caseNum);
    end
end

percent_osc_new_mean = mean(percent_osc_new,1);
conv_in_range_mean = mean(conv_in_range,1);
accuracy1_mean = mean(accuracy1,1);
accuracy2_mean = mean(accuracy2,1);
accuracy3_mean = mean(accuracy3,1);
MSE_testErr_mean = mean(MSE_testErr,1);

percent_osc_new_std = std(percent_osc_new,[],1);
conv_in_range_std = std(conv_in_range,[],1);
accuracy1_std = std(accuracy1,[],1);
accuracy2_std = std(accuracy2,[],1);
accuracy3_std = std(accuracy3,[],1);
MSE_testErr_std = std(MSE_testErr,[],1);

means = [percent_osc_new_mean;...
    conv_in_range_mean;...
    accuracy1_mean;
    accuracy2_mean;
    accuracy3_mean];
stdevs = [percent_osc_new_std;...
    conv_in_range_std;...
    accuracy1_std;
    accuracy2_std;
    accuracy3_std];

% Names = {' ','10neurons',' ','20neurons ',' ','30neurons ',' ','40neurons ',' ','50neurons'};
Names = {' ',' ','10neurons',' ',' ',' ',' ','20neurons ',' '};
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

% % plot %conv and MSE err on a plot with different axes:
% figure;
% xlabel('hidden neuron num');
% yyaxis left
% plot(HiddenN,percent_osc_new_mean);
% ylabel('%inRange');
% 
% yyaxis right
% plot(HiddenN,MSE_testErr_std);
% ylabel('NN MSE error on test group');


% figure;
% plotyy(HiddenN,percent_osc_new_mean,HiddenN,MSE_testErr_std,'bar','plot')

%% Figure 8b:
% non-osc CPGs that produced oscillations after pre-tuning with NN with 'n'
% hidden neurons
% (of different configuration) 

clc 

caseNum = [1,3,5,7,9];
HiddenN = 20;

samplNum = size(sampl,2);
trainRatio = 0.7;
valRatio = 0.15;
testRatio = 1-trainRatio-valRatio;

% find NOT oscilatory CPGs:
load('MatsRandomRes_4Neurons_4Paper_not_osc_Sims.mat','results_not_osc');
howMuch = length(results_not_osc);
cpg_non_osc = randsample(1:howMuch,500);
results_old = results_not_osc(cpg_non_osc);
clear results_not_osc howMuch

seq_n_osc_4Sim = (vertcat(results_old (:).seq))';
seq_n_osc_4Sim = seq_n_osc_4Sim(1:18,:);

numRepeat = 5;

accuracy1 = zeros(numRepeat,length(caseNum));
accuracy2 = zeros(numRepeat,length(caseNum));
accuracy3 = zeros(numRepeat,length(caseNum));
percent_osc_new = zeros(numRepeat,length(caseNum));
conv_in_range = zeros(numRepeat,length(caseNum));
MSE_testErr = zeros(numRepeat,length(caseNum));

for i=1:length(caseNum)
    disp(['NN case #',num2str(caseNum(1,i)),':']);
    
    [inputsNames_4GA,outputsNames_4GA] = ...
    check_NN_case_for_paper(caseNum(1,i),'period_desired');

    [sampl_4GA,targ_4GA] = prepare_NN_inOut(seq_n_osc_4Sim,cpg_non_osc,...
    inputsNames_4GA,outputsNames_4GA,seqOrder);

    [inputsNames_4NN,outputsNames_4NN] = ...
        check_NN_case_for_paper(caseNum(1,i),'period');

    [sampl,targ] = prepare_NN_inOut(seq_in_range,periods_in_range,...
        inputsNames_4NN,outputsNames_4NN,seqOrder);

    
    for j=1:numRepeat
        disp(['    iter ',num2str(j),' out of ',num2str(numRepeat)]);
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
        if 1
            [net,tr] = train(net, sampl, targ);
            MSE_testErr(j,i) = tr.best_tperf;
        else % use untrained net (with random weights)
            net = configure(net,sampl,targ);
            MSE_testErr(j,i) = NaN;
        end
        
        [percent_osc_new(j,i),conv_in_range(j,i),...
            accuracy1(j,i),accuracy2(j,i),accuracy3(j,i)] = ...
             NN_GA_perf(MML,net,sampl_4GA,results_old,seqOrder,caseNum(1,i));
    end
end

percent_osc_new_mean = mean(percent_osc_new,1);
conv_in_range_mean = mean(conv_in_range,1);
accuracy1_mean = mean(accuracy1,1);
accuracy2_mean = mean(accuracy2,1);
accuracy3_mean = mean(accuracy3,1);
MSE_testErr_mean = mean(MSE_testErr,1);

percent_osc_new_std = std(percent_osc_new,[],1);
conv_in_range_std = std(conv_in_range,[],1);
accuracy1_std = std(accuracy1,[],1);
accuracy2_std = std(accuracy2,[],1);
accuracy3_std = std(accuracy3,[],1);
MSE_testErr_std = std(MSE_testErr,[],1);

means = [percent_osc_new_mean;...
    conv_in_range_mean;...
    accuracy1_mean;
    accuracy2_mean;
    accuracy3_mean];
stdevs = [percent_osc_new_std;...
    conv_in_range_std;...
    accuracy1_std;
    accuracy2_std;
    accuracy3_std];

Names = {' ','case#1',' ','case#3',' ','case#5',' ','case#7',' ','case#9'};
label_Y = '';
graph_title = ['NN perf over different case#'];
graph_legend = {'converged','conv in range','accuracy'};
plot_bars_with_errors(means,stdevs,...
    Names,label_Y,graph_title,graph_legend);

disp('|-------------------------------------------------------------|');
disp('| case    | conv mean | conv std | inRange mean | inRange std |');
disp('|-------------------------------------------------------------|');
res = cell(1,length(caseNum));
for i=1:length(caseNum)
    res{1,i} = sprintf('|case %d   | %0.3f     | %0.3f    | %0.3f        | %0.3f       | \n',...
        caseNum(1,i),means(1,i),stdevs(1,i),means(2,i),stdevs(2,i));
    disp(res{1,i});
end
disp('|-------------------------------------------------------------|');

% % plot %conv and MSE err on a plot with different axes:
% figure;
% xlabel('hidden neuron num');
% yyaxis left
% plot(HiddenN,percent_osc_new_mean);
% ylabel('%inRange');
% 
% yyaxis right
% plot(HiddenN,MSE_testErr_std);
% ylabel('NN MSE error on test group');

%% Check avg run time with and without the CB model:
clc

conv_in_range_temp = (periodsMean > MML.perLimOut(1,1)) & ...
    (periodsMean < MML.perLimOut(1,2));

res = results_osc(conv_in_range_temp);
clear conv_in_range_temp

% with:
[avg_sim_time,simOutType] = run_CPG_with_CB(res);

% without:


%% Figure 5: For osc CPGs (in and out of range

clc 

tau_not_in_range = seq_not_in_range(strcmp('tau',seqOrder),:);
tau_in_range = seq_in_range(strcmp('tau',seqOrder),:);
hist_compare(tau_not_in_range,tau_in_range,...
    'tau',20,{'not in range CPGs','in range CPGs'},'plot');
dist_tau = KL_div_4paper(tau_not_in_range,tau_in_range)
clear tau_n_osc tau_osc

b_not_in_range = seq_not_in_range(strcmp('b',seqOrder),:);
b_in_range = seq_in_range(strcmp('b',seqOrder),:);
hist_compare(b_not_in_range,b_in_range,...
    'b',20,{'not in range CPGs','in range CPGs'},'plot');
dist_b = KL_div_4paper(b_not_in_range,b_in_range)
clear b_n_osc b_osc

