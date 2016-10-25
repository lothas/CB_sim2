data1 = load('MatsRandomRes.mat');
data2 = load('MatsScaledRes.mat');
data3 = load('MatsNNRes1.mat');

perRange = [min(data3.desPeriod(:)), max(data3.desPeriod(:))];
randPeriods = data1.periods(~isnan(data1.periods));
randPeriods(randPeriods>perRange(2)*5) = [];
randPeriods(randPeriods<perRange(1)/5) = [];
scalPeriods = data2.periods(~isnan(data2.periods));
scalPeriods(scalPeriods>perRange(2)*5) = [];
scalPeriods(scalPeriods<perRange(1)/5) = [];

% Plot histogram of periods
% figure
% hold on
% histogram(randPeriods,50,'facecolor',[1,0,0],'facealpha',0.5)
% histogram(scalPeriods,50,'facecolor',[0,0,1],'facealpha',0.5)

[N1, bins1] = histcounts(randPeriods,100);
N1 = N1/max(N1);
[N2, bins2] = histcounts(scalPeriods,100);
N2 = N2/max(N2);

bins1n = (bins1 - min(bins1))/(max(bins1) - min(bins1));

% Plot heat map of periods
figure
subplot(2,1,1)
imagesc(N1)
axis off
title('Period distribution of random parameters')
subplot(2,1,2)
imagesc(N2)
set(gca,'YTickLabel',[]);
set(gca,'XTick',bins1n(1:10:end)*100);
set(gca,'XTickLabel',bins1(1:10:end));
title('Period distribution of random parameters after re-scaling')

% Show performance with NN
perf = (data3.netPerf(:,2)+data3.netPerf(:,3))./data3.netPerf(:,4);
bestArc = find(perf == max(perf), 1, 'first');
NNPeriods = data3.sampPerf(bestArc,:);
scNNPeriods = data3.sampPerfSc(bestArc,:);





%%%%%%%%%%%%%%%%%%%%%%%%
netPerf = zeros(1, 4);    % Array to store NN performance
sampPerf = zeros(1, NNSamples);
sampPerfSc = zeros(1, NNSamples);

genomeObj = MML.Gen;
genTauMin = MML.Gen.Range(1,1);
genTauMax = MML.Gen.Range(2,1);
runSimFuncHandle = @MML.runSim;
perLimMin = MML.perLim(1);
perLimMax = MML.perLim(2);

NNSamples = 500;
desPeriod = MML.perLim(1) + ...
             rand(1, NNSamples)*(MML.perLim(2)-MML.perLim(1));
parfor j = 1:NNSamples
    % Setup tau_r, tau_a, c, W and feedback gains using genome
%         seq = getRandFuncHandle(); % Get random genetic sequence
    seq = genomeObj.RandSeq(); %#ok<PFBNS> % Get random genetic sequence

    % Set random b
    if rand()>0.7
        beta = min(max(0.6+0.1*randn(),0.2),0.8);
    else
        beta = min(max(2.5+randn(),0.8),8);
    end

    % Get Tr estimates from NN
    genes = MML.Gen.GetGenes(seq, MML.sample_genes);
    in = [genes';  desPeriod(j)];

    % Normalize X
    in = bsxfun(@rdivide, bsxfun(@minus, in, ...
        MML.normParams(:,1)), MML.normParams(:,2));
    
    seq(1) = predict(Mdl,in');
    if seq(1) < genTauMin || seq(1) > genTauMax
        % Bound tau gene
        seq(1) = min(max(seq(1), genTauMin), genTauMax);
    end

    % Run simulation
    [out, ~, ~] = runSimFuncHandle(seq, beta);    

    % Save resulting period
    if any(isnan(out.periods))
        sampPerf(j) = NaN;
    else
        sampPerf(j) = max(out.periods);
    end

    % Did period converge to the desired range?
    if sampPerf(j) < perLimMin || sampPerf(j) > perLimMax
        % Try rescaling time
        ratio = desPeriod(j)/sampPerf(j);
        seq(1) = seq(1)*ratio;
        if seq(1) < genTauMin || seq(1) > genTauMax
            % Bound tau gene
            seq(1) = min(max(seq(1), genTauMin), genTauMax);
        end

        % Run simulation
        [out, ~, ~] = runSimFuncHandle(seq, beta);
        % Save resulting period
        if any(isnan(out.periods))
            sampPerfSc(j) = NaN;
        else
            sampPerfSc(j) = max(out.periods);
        end
    else
        sampPerfSc(j) = sampPerf(j);
    end
end
