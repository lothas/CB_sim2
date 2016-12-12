MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.01; % 0.01
MML.tEnd = 15; % 15

load('MatsRandomRes.mat');

% Get genomes that converged
ids_conv = find(~isnan(periods));
disp(['Converged: ', int2str(length(ids_conv)), ' out of ', ...
    int2str(length(periods))])

% Get genomes that converged and have a period up to 5 seconds
ids_conv = find(~isnan(periods) & periods<5);
disp(['Converged: ', int2str(length(ids_conv)), ' out of ', ...
    int2str(length(periods))])

sPeriods = periods(ids_conv);
sSeqs = vertcat(results(ids_conv).seq);

%% Permutations
% Chose fewer samples at random
nSamples = 10000;
ids = randsample(1:length(sPeriods),nSamples);
fPeriods = sPeriods(ids);
fSeqs = sSeqs(ids,:);

% Create permutations
permSeqs = MML.permSeq(fSeqs(1,:));
nPerms = size(permSeqs,1);
permSeqs = zeros(nPerms*size(seqs,1), size(seqs,2));
for i = 1:nSamples
    permSeqs((i-1)*nPerms+1:i*nPerms,:) = MML.permSeq(fSeqs(i,:));
end
permPers = reshape(repmat(fPeriods,nPerms,1),[],1);

% Prepare samples and targets
samples = permSeqs(:,[1,2,7:18]);
targets = permPers;

%% Symmetric functions

targets = sPeriods';
samples = zeros(length(targets),9);
eps = 1e-4;
one_over = 0;

for s = 1:length(targets)
    weights = sSeqs(s,7:18);
    
    W = zeros(4, 4);
    k = 1;
    for i = 1:4
        for j = 1:4
            if i == j
                continue
            end
            W(i,j) = weights(k);
            k = k+1;
        end
    end
    
    res1 = 0;
    res2 = 1;
    for m = 1:4
        for n = m+1:4
            res1 = res1 + W(m,n)*W(n,m);
            res2 = res2 * (W(m,n)+W(n,m));
        end
    end
    
    samples(s,1) = sSeqs(s,1); % tau
    if one_over == 1
        samples(s,2) = 1/sSeqs(s,2); % b
        % srqt(Product of all weights)
        samples(s,3) = 1/(eps + sqrt(abs(prod(weights))));
        % Sum of all weights
        samples(s,4) = 1/sum(weights);
        % Prod of sum of weight going into each neuron (close to 0)
        samples(s,5) = 1/(eps + prod(sum(reshape(weights,3,4))));
        % Sum of prod of weight going into each neuron (close to 0)
        samples(s,6) = 1/(eps + sum(prod(reshape(weights,3,4),1)));
        % Sum of prod of weight pairs (close to 0)
        samples(s,7) = 1/(eps + res1);
        % Prod of sum of weight pairs (close to 0)
        samples(s,8) = 1/(eps + res2);
    else
        samples(s,2) = sSeqs(s,2); % b
        % srqt(Product of all weights)
        samples(s,3) = sqrt(abs(prod(weights)));
        % Sum of all weights
        samples(s,4) = sum(weights);
        % Prod of sum of weight going into each neuron (close to 0)
        samples(s,5) = prod(sum(reshape(weights,3,4)));
        % Sum of prod of weight going into each neuron (close to 0)
        samples(s,6) = sum(prod(reshape(weights,3,4),1));
        % Sum of prod of weight pairs (close to 0)
        samples(s,7) = res1;
        % Prod of sum of weight pairs (close to 0)
        samples(s,8) = res2;
    end
    % Determinant of weight matrix
    samples(s,9) = det(W);
end

% Remove "outliers"
sampMean = mean(samples);
sampStd = std(samples);
goodIds = ones(size(samples,1),1);
for i = 1:9
    goodIds = goodIds & samples(:,i) < sampMean(i)+2*sampStd(i) ...
        & samples(:,i) > sampMean(i)-2*sampStd(i);
end
gSamples = samples(goodIds,:);
gTargets = targets(goodIds,:);

figure
for p = 1:9
    subplot(3,3,p)
    hist(samples(:,p),100);
end
        
% Normalize data
gsampMean = mean(gSamples);
gsampStd = std(gSamples);
nSamples = bsxfun(@rdivide,bsxfun(@minus, gSamples, gsampMean), gsampStd);

figure
for p = 1:9
    subplot(3,3,p)
    hist(nSamples(:,p),100);
end