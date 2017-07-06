%% Initialize machine learning object for Matsuoka analysis
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01; % 0.01
MML.tEnd = 15; % 15

% % Set constant beta
% MML.Sim.Con.beta = 7;

nPlotSamples = 0; % 10;

% Turn off findpeaks warning
warning('off','signal:findpeaks:largeMinPeakHeight');

%% Test effect of permutating genes
while 1
    seq = MML.Gen.RandSeq(1);
    [out, ~, ~] = MML.runSim(seq);
    if all(~isnan(out.periods))
        break
    end
end
permSeqs = MML.permSeq(seq);
nPerms = size(permSeqs,1);
out = cell(1,nPerms); sim = cell(1,nPerms); signal = cell(1,nPerms);
for i = 1:nPerms
    [out{i}, sim{i}, signal{i}] = MML.runSim(permSeqs(i,:));
end

% Show results
figure
hold on
periods = zeros(nPerms,1);
perSignal = signal;
for i = 1:nPerms
    if any(isnan(out{i}.periods))
        periods(i) = NaN;
    else
        periods(i) = max(out{i}.periods);
    end
    
    % Get last N cycles
    [sigT, perSignal{i}] = MML.getNcycles(signal{i}.T, ...
        signal{i}.signal, periods(i), 5);
    plot(sigT, perSignal{i});
end
disp(['Converged: ', int2str(sum(~isnan(periods))), ...
    ' out of ', int2str(nPerms)])
    

%% Phase 1 - Run lots of Matsuoka simulations with different parameters
filename1 = 'MatsRandomRes.mat';
% filename1 = 'MatsRandomRes_test.mat';
nSamples = 200000;
MML.runRandomSims(nSamples, filename1);

%% Phase 1.1 - Show period calculation errors

data = load(filename1);
conv = true(data.nSims,1);
for i = 1:data.nSims
    if any(isnan(data.results(i).periods))
        conv(i) = false;
    end
end
id_conv = find(conv == true);
id_not_conv = find(conv == false);
disp(['Detected period in ', int2str(length(id_conv)), ' out of ', ...
    int2str(data.nSims), ' simulations (', ...
    num2str(length(id_conv)/data.nSims*100,'%.1f'), '%)'])

perOK1 = vertcat(data.results.perOK1);
perError1_1 = zeros(length(id_conv), 2*MML.nNeurons);
nw = 0;
for i = 1:length(id_conv)
    if length(data.results(id_conv(i)).perError1) ~= 2*MML.nNeurons
        disp(id_conv(i))
        nw = nw+1;
    end
    perError1_1(i,:) = data.results(id_conv(i)).perError1;
end
perError1_2 = horzcat(data.results(id_not_conv).perError1)';

perOK2 = vertcat(data.results.perOK2);
perError2_1 = horzcat(data.results(id_conv).perError2)';
perError2_2 = horzcat(data.results(id_not_conv).perError2)';
% When no period was detected Error2 is always NaN.
% So perError2_2 is useless, but perOK2 still contains helpful info
% (the threshold used was also 1e-6)

perError1_1_thresh = 0.1;
disp([num2str(sum(any(perError1_1'<perError1_1_thresh)) ...
    / length(id_conv) * 100, '%.1f'), ...
    '% of oscillatory CPGS have period detection error under ', ...
    num2str(perError1_1_thresh)])

perError1_2_thresh = 1e-6;
disp([num2str(sum(any(perError1_2'<perError1_2_thresh)) ...
    / length(id_not_conv) * 100, '%.1f'), ...
    '% of non oscillatory CPGS have period detection error under ', ...
    num2str(perError1_2_thresh)])

perError2_1_thresh = 0.0001;
disp([num2str(sum(any(perError2_1'<perError2_1_thresh)) ...
    / length(id_conv) * 100, '%.1f'), ...
    '% of oscillatory CPGS have period detection error under ', ...
    num2str(perError2_1_thresh)])

% show histograms
perError1_1(any(perError1_1'>2*perError1_1_thresh),:) = [];
perError1_2(any(perError1_2'>2*perError1_2_thresh),:) = [];
perError2_1(any(perError2_1'>2*perError2_1_thresh),:) = [];
figure
subplot(1,3,1)
hist(perError1_1(:),50)
title('Error of oscillatory')
subplot(1,3,2)
hist(perError1_2(:),50)
title('Error of non oscillatory')
subplot(1,3,3)
hist(perError2_1(:),50)
title('Error of oscillatory (Rea)')

neuronAct = vertcat(data.results.neuronActive);
neuronOsc = vertcat(data.results.neuronOsc);

% Definitely wrong:
% Only one neuron is oscillatory
ids1 = find(sum(neuronOsc,2)==1);
% One output signal is 0 (pair of neurons not oscillatory)
ids2 = find(neuronOsc(:,1)==0 & neuronOsc(:,3)==0 | neuronOsc(:,2)==0 & neuronOsc(:,4)==0);

% Could be OK:
ids3 = find(neuronOsc(:,1)==1 & neuronOsc(:,2)==1 | neuronOsc(:,3)==1 & neuronOsc(:,4)==1);

% Problematic detections (different period for out1 and out2)
ids4 = find(all(~isnan(per')) & diff(per')>1e-3);

% Sum it all up
detOsc = false(data.nSims,1); % Oscillations detected
perOK1 = false(data.nSims,1); % Period calculation error 1
perOK2 = false(data.nSims,1); % Period calculation error 2
perOK3 = false(data.nSims,1); % Periodic neurons
for i = 1:data.nSims
    detOsc(i) = all(~isnan(data.results(i).periods));
    if detOsc(i)
        perOK1(i) = max(data.results(i).perError1)<perError1_1_thresh;
        perOK2(i) = max(data.results(i).perError2)<perError2_1_thresh;
    else
        % No period was detected
        perOK1(i) = max(data.results(i).perError1)<perError1_2_thresh;
        perOK2(i) = all(data.results(i).perOK2);
    end
    
    oscN = data.results(i).neuronOsc;
    perOK3(i) = true;
    % Good combinations
    for j = 1:MML.nNeurons/2
        % If these pairs are non-oscillatory, their corresponding output
        % will be stationary
        perOK3(i) = perOK3(i) & (oscN(2*(j-1)+1) | oscN(2*j));
    end
end

disp(['Detected period in ', int2str(sum(detOsc)), ' out of ', ...
    int2str(data.nSims), ' simulations (', ...
    num2str(sum(detOsc)/data.nSims*100,'%.1f'), '%)'])
disp(['Should have period based on detected oscillatory neurons: ', ...
    int2str(sum(perOK3)), ' out of ', int2str(data.nSims), ' simulations (', ...
    num2str(sum(perOK3)/data.nSims*100,'%.1f'), '%)'])
disp(['Correctly detected: ', int2str(sum(detOsc & perOK3))])
disp(['Falsely detected: ', int2str(sum(detOsc & ~perOK3))])
disp(['Missed: ', int2str(sum(~detOsc & perOK3))])
disp(' ')

disp('Correctly detected oscillations:')
disp(['* Jonathan method: ', int2str(sum(perOK1)), ' out of ', ...
   int2str(data.nSims), ' simulations (', ...
    num2str(sum(perOK1)/data.nSims*100,'%.1f'), '%)'])
disp(['* Rea method: ', int2str(sum(perOK2)), ' out of ', ...
   int2str(data.nSims), ' simulations (', ...
    num2str(sum(perOK2)/data.nSims*100,'%.1f'), '%)'])
disp(['* Both methods: ', int2str(sum(perOK1 & perOK2)), ' out of ', ...
   int2str(data.nSims), ' simulations (', ...
    num2str(sum(perOK1 & perOK2)/data.nSims*100,'%.1f'), '%)'])
disp(' ')

% We define good results as those confirmed by both methods,
% correctly detected (period detected and neurons are oscillatory),
% and all signals have the same period
periods = horzcat(data.results.periods);
good_ids = find(perOK1 & perOK2 & ...
    (detOsc & perOK3 | ~detOsc & ~perOK3) & ...
    (~detOsc | abs(diff(periods))' < 0.5*min(periods(:))));
results = data.results(good_ids);
nSims = length(results);
gOsc = detOsc(good_ids);
gPer = perOK3(good_ids);

disp(['Detected period in ', int2str(sum(gOsc)), ' out of ', ...
    int2str(nSims), ' simulations (', ...
    num2str(sum(gOsc)/nSims*100,'%.1f'), '%)'])
disp(['Should have period based on detected oscillatory neurons: ', ...
    int2str(sum(gPer)), ' out of ', int2str(nSims), ' simulations (', ...
    num2str(sum(gPer)/nSims*100,'%.1f'), '%)'])
disp(['Correctly detected: ', int2str(sum(gOsc & gPer))])
disp(['Falsely detected: ', int2str(sum(gOsc & ~gPer))])
disp(['Missed: ', int2str(sum(~gOsc & gPer))])
disp(' ')

periods = horzcat(results.periods)';
% speriods = sortrows(periods,-1);
% plot(speriods(:,1),speriods(:,2))
periods = max(periods,[],2);
id_conv = find(gOsc);
id_per = find(periods >= MML.perLimOut(1) & periods <= MML.perLimOut(2));
disp(['Produced desired period range: ', int2str(length(id_per)), ' out of ', ...
    int2str(nSims), ' simulations (', ...
    num2str(length(id_per)/nSims*100,'%.1f'), '%)'])
save('MatsRandomChkRes.mat','nSims','results','id_conv','id_per','periods');
data.nSims = nSims;
data.results = results;
data.periods = max(periods,[],2);

%% Phase 2 - Re-run simulations that converged outside the desired range,
% this time with scaled temporal parameters
filename2 = 'MatsScaledRes.mat';
% data = load(filename1);
reDo_ids = zeros(1, data.nSims);
reDo_ids(~isnan(data.periods)) = 1;
reDo_ids(data.periods >= MML.perLimOut(1) & ...
    data.periods <= MML.perLimOut(2)) = 0;
% reDo_ids(data.id_conv) = 1;
% reDo_ids(data.id_per) = 0;
inputData = data.results(logical(reDo_ids));
inputPeriods = data.periods(logical(reDo_ids));
MML.runScaledSims(inputData, inputPeriods, filename2);

%% Phase 2.1 - check results

data = load(filename2);
detOsc = false(data.nSims,1); % Oscillations detected
perOK1 = false(data.nSims,1); % Period calculation error 1
perOK2 = false(data.nSims,1); % Period calculation error 2
perOK3 = false(data.nSims,1); % Periodic neurons
for i = 1:data.nSims
    detOsc(i) = all(~isnan(data.results(i).periods));
    if detOsc(i)
        perOK1(i) = max(data.results(i).perError1)<perError1_1_thresh;
        perOK2(i) = max(data.results(i).perError2)<10*perError2_1_thresh;
    else
        % No period was detected
        perOK1(i) = max(data.results(i).perError1)<perError1_2_thresh;
        perOK2(i) = all(data.results(i).perOK2);
    end
    
    oscN = data.results(i).neuronOsc;
    perOK3(i) = true;
    % Good combinations
    for j = 1:MML.nNeurons/2
        % If these pairs are non-oscillatory, their corresponding output
        % will be stationary
        perOK3(i) = perOK3(i) & (oscN(2*(j-1)+1) | oscN(2*j));
    end
end

disp(['Detected period in ', int2str(sum(detOsc)), ' out of ', ...
    int2str(data.nSims), ' simulations (', ...
    num2str(sum(detOsc)/data.nSims*100,'%.1f'), '%)'])
disp(['Should have period based on detected oscillatory neurons: ', ...
    int2str(sum(perOK3)), ' out of ', int2str(data.nSims), ' simulations (', ...
    num2str(sum(perOK3)/data.nSims*100,'%.1f'), '%)'])
disp(['Correctly detected: ', int2str(sum(detOsc & perOK3))])
disp(['Falsely detected: ', int2str(sum(detOsc & ~perOK3))])
disp(['Missed: ', int2str(sum(~detOsc & perOK3))])
disp(' ')

disp('Correctly detected oscillations:')
disp(['* Jonathan method: ', int2str(sum(perOK1)), ' out of ', ...
   int2str(data.nSims), ' simulations (', ...
    num2str(sum(perOK1)/data.nSims*100,'%.1f'), '%)'])
disp(['* Rea method: ', int2str(sum(perOK2)), ' out of ', ...
   int2str(data.nSims), ' simulations (', ...
    num2str(sum(perOK2)/data.nSims*100,'%.1f'), '%)'])
disp(['* Both methods: ', int2str(sum(perOK1 & perOK2)), ' out of ', ...
   int2str(data.nSims), ' simulations (', ...
    num2str(sum(perOK1 & perOK2)/data.nSims*100,'%.1f'), '%)'])
disp(' ')

% We define good results as those confirmed by both methods,
% correctly detected (period detected and neurons are oscillatory),
% and all signals have the same period
periods = horzcat(data.results.periods);
good_ids = find(perOK1 & perOK2 & ...
    (detOsc & perOK3 | ~detOsc & ~perOK3) & ...
    (~detOsc | abs(diff(periods))' < 0.1*min(periods(:))));
results = data.results(good_ids);
nSims = length(results);
gOsc = detOsc(good_ids);
gPer = perOK3(good_ids);

disp(['Detected period in ', int2str(sum(gOsc)), ' out of ', ...
    int2str(nSims), ' simulations (', ...
    num2str(sum(gOsc)/nSims*100,'%.1f'), '%)'])
disp(['Should have period based on detected oscillatory neurons: ', ...
    int2str(sum(gPer)), ' out of ', int2str(nSims), ' simulations (', ...
    num2str(sum(gPer)/nSims*100,'%.1f'), '%)'])
disp(['Correctly detected: ', int2str(sum(gOsc & gPer))])
disp(['Falsely detected: ', int2str(sum(gOsc & ~gPer))])
disp(['Missed: ', int2str(sum(~gOsc & gPer))])
disp(' ')

periods = horzcat(results.periods)';
speriods = sortrows(periods,-1);
plot(speriods(:,1),speriods(:,2))
periods = max(periods,[],2);
id_conv = find(gOsc);
id_per = find(periods >= MML.perLimOut(1) & periods <= MML.perLimOut(2));
disp(['Produced desired period range: ', int2str(length(id_per)), ' out of ', ...
    int2str(nSims), ' simulations (', ...
    num2str(length(id_per)/nSims*100,'%.1f'), '%)'])
save('MatsScaledChkRes.mat','nSims','results','id_conv','id_per','periods');
data.nSims = nSims;
data.results = results;
data.periods = max(periods,[],2);

%% Phase 2.2 - Check damping condition
% data = load(filename2);
results = data.results;
converged = ~isnan(data.periods);
passed_cond1 = zeros(length(results),1);
tp1 = passed_cond1; fp1 = tp1; tn1 = fp1; fn1 = tp1;
for i = 1:length(results)
    cr = results(i);
%     disp(['(',num2str(cr.Tr,5),' - ',num2str(cr.Ta,5),')^2 = ',...
%         num2str((cr.Tr-cr.Ta)^2,5)]);
    if (cr.Tr-cr.Ta)^2 >= 4*cr.Tr*cr.Ta*cr.b
        passed_cond1(i) = 1;
%         disp('              >=');
    else
%         disp('              <');
    end
%     disp(['4*',num2str(cr.Tr,5),'*',num2str(cr.Ta,5),'*',num2str(cr.b,5),' = ',...
%         num2str(4*cr.Tr*cr.Ta*cr.b,5)]);
%     disp(' ');
    
    if passed_cond1(i)
        if converged(i)
            tp1(i) = 1;
        else
            fp1(i) = 1;
        end
    else
        if converged(i)
            fn1(i) = 1;
        else
            tn1(i) = 1;
        end
    end
end
% Show results
disp(['Converged: ', int2str(sum(converged)), ...
    ' out of ', int2str(numel(converged))]);
disp(['Passed condition 1: ', int2str(sum(passed_cond1)), ...
    ' out of ', int2str(numel(passed_cond1))]);
disp(['True positives: ', int2str(sum(tp1))]);
disp(['True negatives: ', int2str(sum(tn1))]);
disp(['False positives: ', int2str(sum(fp1))]);
disp(['False negatives: ', int2str(sum(fn1))]);

% Show some cases
n_cases = nPlotSamples;
% True positives
MML.plotSamples(results, tp1, n_cases, 'Cond. 1 true positive sample #');
% True negatives
MML.plotSamples(results, tn1, n_cases, 'Cond. 1 true negative sample #');
% False positives
MML.plotSamples(results, fp1, n_cases, 'Cond. 1 false positive sample #');
% False negatives
MML.plotSamples(results, fn1, n_cases, 'Cond. 1 false negative sample #');

%% Phase 2.2 - Check tonic input condition
passed_cond2 = zeros(length(results),length(results(1).c));
tp2 = 0*converged; fp2 = tp2; tn2 = fp2; fn2 = tp2;
for i = 1:length(results)
    cr = results(i);
    passed_cond2(i, :) = (cr.c > cr.W*cr.c/(1+cr.b))';
        
    if all(passed_cond2(i, :))
        if converged(i)
            tp2(i) = 1;
        else
            fp2(i) = 1;
        end
    else
        if converged(i)
            fn2(i) = 1;
        else
            tn2(i) = 1;
        end
    end
end
% Show results
disp(['Converged: ', int2str(sum(converged)), ...
    ' out of ', int2str(numel(converged))]);
disp(['Passed condition 2: ', int2str(sum(all(passed_cond2'))), ...
    ' out of ', int2str(numel(converged))]);
disp(['True positives: ', int2str(sum(tp2))]);
disp(['True negatives: ', int2str(sum(tn2))]);
disp(['False positives: ', int2str(sum(fp2))]);
disp(['False negatives: ', int2str(sum(fn2))]);

% Show some cases
n_cases = nPlotSamples;
% True positives
MML.plotSamples(results, tp2, n_cases, 'Cond. 2 true positive sample #');
% True negatives
MML.plotSamples(results, tn2, n_cases, 'Cond. 2 true negative sample #');
% False positives
MML.plotSamples(results, fp2, n_cases, 'Cond. 2 false positive sample #');
% False negatives
MML.plotSamples(results, fn2, n_cases, 'Cond. 2 false negative sample #');

%% Phase pre-3 - Train NNs using different features and targets

% Combinations for genome with tau, beta, amp and weights.
sample_genes = {{'amp','weights'};
                {'beta','amp','weights'};
                {'weights'};
                {'beta','weights'};
                {'\tau_r','weights'};
                {'\tau_r','amp','weights'}};
target_genes = {{'\tau_r'};
                {'\tau_r','beta'};
                {'beta'}};
combos = [1,1; % {'amp','weights'}          -> {'\tau_r'}           0.52    0.18
          1,2; % {'amp','weights'}          -> {'\tau_r','beta'}    0.81    0.14 ++
          1,3; % {'amp','weights'}          -> {'beta'}             0.83    0.10 +gr
          2,1; % {'beta','amp','weights'}   -> {'\tau_r'}           0.52    0.18
          3,1; % {'weights'}                -> {'\tau_r'}           0.54    0.20
          3,2; % {'weights'}                -> {'\tau_r','beta'}    0.78    0.12 +++
          3,3; % {'weights'}                -> {'beta'}             0.80    0.06 +
          4,1; % {'beta','weights'}         -> {'\tau_r'}           0.50    0.20 
          5,3; % {'\tau_r','weights'}       -> {'beta'}             0.85    0.09 +
          6,3];% {'\tau_r','amp','weights'} -> {'beta'}             0.85    0.09 +
      
% Combinations for genome with tau, tau_ratio, beta, amp and weights.
% sample_genes = {{'amp','weights'};
%                 {'beta','amp','weights'};
%                 {'weights'};
%                 {'beta','weights'};
%                 {'\tau_r','weights'};
%                 {'amp','weights','\tau_ratio'};
%                 {'beta','amp','weights','\tau_ratio'};
%                 {'weights','\tau_ratio'};
%                 {'beta','weights','\tau_ratio'};
%                 {'\tau_r','weights','\tau_ratio'}};
% target_genes = {{'\tau_r'};
%                 {'\tau_r','beta'};
%                 {'beta'};
%                 {'\tau_r','\tau_ratio'};
%                 {'\tau_r','beta','\tau_ratio'};
%                 {'beta','\tau_ratio'}};
% combos = [1,1; % {'amp','weights'}          -> {'\tau_r'}           0.50    0.19
%           1,2; % {'amp','weights'}          -> {'\tau_r','beta'}    0.80    0.14 ++
%           1,3; % {'amp','weights'}          -> {'beta'}             0.81    0.08 +
%           1,4;
%           1,5;
%           1,6; % **** 6
%           2,1; % {'beta','amp','weights'}   -> {'\tau_r'}           0.49    0.16
%           2,4;
%           3,1; % {'weights'}                -> {'\tau_r'}           0.52    0.18
%           3,2; % **** 10 {'weights'}                -> {'\tau_r','beta'}    0.81    0.18 +++
%           3,3; % **** 11 {'weights'}                -> {'beta'}             0.79    0.09 +
%           3,4;
%           3,5; % **** 13
%           3,6; % **** 14
%           4,1; % {'beta','weights'}         -> {'\tau_r'}           0.52    0.19 
%           4,4;
%           5,3; % **** 17 {'\tau_r','weights'}      -> {'beta'}             0.83    0.08 +
%           5,6;
%           6,1; % {'amp','weights'}          -> {'\tau_r'}           0.50    0.19
%           6,2; % **** 20 {'amp','weights'}          -> {'\tau_r','beta'}    0.80    0.14 ++
%           6,3; % **** 21{'amp','weights'}          -> {'beta'}             0.81    0.08 +
%           7,1; % {'beta','amp','weights'}   -> {'\tau_r'}           0.49    0.16
%           8,1; % {'weights'}                -> {'\tau_r'}           0.52    0.18
%           8,2; % {'weights'}                -> {'\tau_r','beta'}    0.81    0.18 +++
%           8,3; % **** 25 {'weights'}                -> {'beta'}             0.79    0.09 +
%           9,1; % {'beta','weights'}         -> {'\tau_r'}           0.52    0.19 
%           10,3]; % **** 27 {'\tau_r','weights'}      -> {'beta'}             0.83    0.08 +
      
% Combinations for genome with tau, amp and weights (no beta)
% sample_genes = {{'amp','weights'};
%                 {'weights'}};
% target_genes = {{'\tau_r'};
%                 {'\tau_r','amp'};
%                 {'amp'}};
% combos = [1,1; % {'amp','weights'}          -> {'\tau_r'}           0.50    0.19
%           2,1; % {'weights'}                -> {'\tau_r'}           0.52    0.18
%           2,2; % {'weights'}                -> {'\tau_r','amp'}    0.81    0.18 +++
%           2,3]; % {'weights'}                -> {'amp'}             0.79    0.09 +

maxN = min(nSamples, 250000);
NNSamples = 500;
% inFilenames = {filename1, filename2};
inFilenames = {'MatsRandomChkRes.mat', 'MatsScaledChkRes.mat'};

nCombos = size(combos,1);
net = cell(nCombos, 1);           % Cell array to store NNs
tr = cell(nCombos, 1);            % Cell array to store NNs training res
netPerf = zeros(nCombos, 4);      % Array to store NN performance
% Array to store NN per sample performance
desPeriod = zeros(nCombos, NNSamples);
sampPerf = zeros(nCombos, NNSamples);
sampPerfSc = zeros(nCombos, NNSamples); % re-scaled results
    
for i = 1:nCombos
    MML.sample_genes = sample_genes{combos(i,1)};
    MML.target_genes = target_genes{combos(i,2)};

    % Normalize weights
%     weight_id = [];
    weight_id = find(strcmp(MML.sample_genes,'weights'),1,'first');
    MML.norm_weights = 0;
    if ~isempty(weight_id)
        MML.norm_weights = 1;
        MML.sample_genes1 = {'weights'};
        MML.sample_genes2 = MML.sample_genes;
        MML.sample_genes2(weight_id) = [];
    end
    
    MML.normParams = [];
    [samples, targets, normParams] = MML.prepareNNData(inFilenames, maxN);
    MML.normParams = normParams;
    
    [net{i}, tr{i}, netPerf(i,:), desPeriod(i,:), ...
            sampPerf(i,:), sampPerfSc(i,:)] = ...
            MML.trainNN(samples, targets, 50, NNSamples);
        
%     [net1, tr1, netPerf1, desPeriod1, sampPerf1, sampPerfSc1] = ...
%             MML.trainNN(samples, targets, 20, NNSamples);
end

save('MatsNNTests2', 'sample_genes', 'target_genes', 'combos', ...
     'nCombos', 'net', 'tr', 'netPerf', 'desPeriod', ...
     'sampPerf', 'sampPerfSc');
 
 % Display results
 for i = 1:nCombos
     disp([num2cell(netPerf(i,:)), ...
         'Feat:', cellstr(sample_genes{combos(i,1)}), ...
         'Target:', cellstr(target_genes{combos(i,2)})])
 end
 
%% Phase 3 - Train NNs using the data from phases 1 and 2

% Pick best features
score = netPerf(:,2) + netPerf(:,3);
best_comb = find(score == max(score));
if length(best_comb) > 1
    best_scores = netPerf(best_comb,2);
    id = find(best_scores == max(best_scores), 1, 'first');
    best_best = best_comb(id);
    best_comb = best_best;
end
disp('Best features/target combo:');
disp(['Features', sample_genes{combos(best_comb,1)}]);
disp(['Targets', target_genes{combos(best_comb,2)}]);
disp(netPerf(best_comb,:));

MML.sample_genes = sample_genes{combos(best_comb,1)};
MML.target_genes = target_genes{combos(best_comb,2)};
    
filename3 = 'MatsNNData.mat';
filename4 = 'MatsNNRes1.mat';
filename5 = 'MatsNNRes2.mat';
filename6 = 'MatsNNRes3.mat';
if exist(filename3, 'file') ~= 2
    maxN = min(nSamples, 250000);
    NNSamples = 500;
    inFilenames = {filename1, filename2};
    [samples, targets, normParams] = MML.prepareNNData(inFilenames, maxN);
    save(filename3, 'maxN', 'NNSamples', ...
                    'samples', 'targets', 'normParams');
else
    load(filename3);
end
MML.normParams = normParams;

% Train networks with different architectures and all samples
if exist(filename4, 'file') ~= 2
    architectures = {5, 20, 50, [25, 25], [15, 25, 10]};
    nArch = numel(architectures);
    
    net = cell(nArch, 1);           % Cell array to store NNs
    tr = cell(nArch, 1);            % Cell array to store NNs training res
    netPerf = zeros(nArch, 4);      % Array to store NN performance
    % Array to store NN per sample performance
    desPeriod = zeros(nArch, NNSamples);
    sampPerf = zeros(nArch, NNSamples);
    sampPerfSc = zeros(nArch, NNSamples); % re-scaled results
    
    for i = 1:nArch
        [net{i}, tr{i}, netPerf(i,:), desPeriod(i,:), ...
            sampPerf(i,:), sampPerfSc(i,:)] = ...
            MML.trainNN(samples, targets, architectures{i}, NNSamples);
    end
    save(filename4, 'architectures', 'nArch', 'net', 'tr', ...
                    'netPerf', 'desPeriod', 'sampPerf', 'sampPerfSc');
else
    load(filename4);
end

% Train network with specific architecture and different number of samples
if exist(filename5, 'file') ~= 2
    maxN = size(samples, 2);
    architecture = architectures{2};
    nSampleSizes = 10;
    sampleSizes = floor(logspace(1,3,nSampleSizes)*maxN/1000);
    
    net = cell(nSampleSizes, 1);       % Cell array to store NNs
    tr = cell(nArch, 1);            % Cell array to store NNs training res
    netPerf = zeros(nSampleSizes, 4);  % Array to store NN performance
    % Array to store NN per sample performance
    desPeriod = zeros(nSampleSizes, NNSamples);
    sampPerf = zeros(nSampleSizes, NNSamples);
    sampPerfSc = zeros(nArch, NNSamples); % re-scaled results
    
    for i = 1:nSampleSizes
        ids = randsample(maxN, sampleSizes(i));
        [net{i}, tr{i}, netPerf(i,:), desPeriod(i,:), ...
            sampPerf(i,:), sampPerfSc(i,:)] = ...
            MML.trainNN(samples(:, ids), targets(:, ids), ...
                        architecture, NNSamples);
    end
    save(filename5, 'nSampleSizes', 'sampleSizes', 'architecture',...
                    'net', 'tr', 'netPerf', 'desPeriod', ...
                    'sampPerf', 'sampPerfSc');
else
    load(filename5);
end

% Train network with specific architecture and different number of samples
if exist(filename6, 'file') ~= 2
    maxN = size(samples, 2);
    architecture = architectures{end};
    nSampleSizes = 10;
    sampleSizes = floor(logspace(1,3,nSampleSizes)*maxN/1000);
    
    net = cell(nSampleSizes, 1);       % Cell array to store NNs
    tr = cell(nArch, 1);            % Cell array to store NNs training res
    netPerf = zeros(nSampleSizes, 4);  % Array to store NN performance
    % Array to store NN per sample performance
    desPeriod = zeros(nSampleSizes, NNSamples);
    sampPerf = zeros(nSampleSizes, NNSamples);
    sampPerfSc = zeros(nArch, NNSamples); % re-scaled results
    
    for i = 1:nSampleSizes
        ids = randsample(maxN, sampleSizes(i));
        [net{i}, tr{i}, netPerf(i,:), desPeriod(i,:), ...
            sampPerf(i,:), sampPerfSc(i,:)] = ...
            MML.trainNN(samples(:, ids), targets(:, ids), ...
                        architecture, NNSamples);
    end
    save(filename6, 'nSampleSizes', 'sampleSizes', 'architecture',...
                    'net', 'tr', 'netPerf', 'desPeriod', ...
                    'sampPerf', 'sampPerfSc');
else
    load(filename6);
end

%% Phase 4 - Train SVMs using the data from phases 1 and 2

%% Phase 5 - ???

%% Phase 6 - Profit! (Plot results)
data1 = load(filename1); % Random results
data2 = load(filename2); % Scaled results
data3 = load(filename4); % NN results - diff. architectures
data4 = load(filename5); % NN results - diff. number of training samples
data5 = load(filename6); % NN results - diff. number of training samples
load(filename3);

MML.plotNNConv(data1, data2, data3, 1);
MML.plotNNConv(data1, data2, data4, 2);
MML.plotNNConv(data1, data2, data5, 2);

nBins = 50;
nDists = 1+size(data3.sampPerf,1);
rows = max(floor(sqrt(nDists)),1);
cols = ceil(nDists/rows);

% figure
% % Plot distribution of random samples
% subplot(rows, cols, 1);
% hist(max([data1.results(data1.id_conv).periods]),nBins)
% for i = 2:nDists
%     % Plot distribution of NN samples
%     subplot(rows, cols, i)
%     hist(data3.sampPerf(i-1,:),nBins)
% end

figure
% Plot distribution of random samples
subplot(rows, cols, 1);
% hist(max([data1.results(data1.id_conv).periods]),nBins)
[counts,centers]=hist(max([data1.results(data1.id_conv).periods]),nBins);
% b1=bar(centers,counts/trapz(centers,counts));
b1=bar(centers,counts/max(counts));
hold on
[counts,centers]=hist(max([data2.results(data2.id_conv).periods]),nBins);
% b2=bar(centers,counts/trapz(centers,counts));
b2=bar(centers,counts/max(counts));
set(get(b1,'Children'),'Facecolor',[0 0 1],'EdgeColor','k','FaceAlpha',0.5);
set(get(b2,'Children'),'Facecolor',[1 0 0],'EdgeColor','k','FaceAlpha',0.5);
% h = findobj(gca,'Type','patch');
% set(h(1),'Facecolor',[1 0 0],'EdgeColor','k','FaceAlpha',0.5);
% set(h(2),'Facecolor',[0 0 1],'EdgeColor','k','FaceAlpha',0.5);

for i = 2:nDists
    % Plot distribution of NN samples
    subplot(rows, cols, i)
    hist(data3.sampPerf(i-1,:),nBins)
end
