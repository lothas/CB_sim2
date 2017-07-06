function [ samples, targets, normParams ] = ...
    prepareNNData(obj, filenames, maxN)
%PREPARENNDATA Prepares data to train a neural network using stored results
    % Load data from files
    data = load(filenames{1});
    
    % Keep only results that converged
    new_periods = data.periods;
    results = data.results(~isnan(new_periods));
    periods = new_periods(~isnan(new_periods));
    for i = 2:numel(filenames)
        data = load(filenames{i});
        new_periods = data.periods;
        results = [results, data.results(~isnan(new_periods))]; %#ok<AGROW>
        periods = [periods, new_periods(~isnan(new_periods))]; %#ok<AGROW>
    end
    nSamples = numel(results);
    
    if nargin < 3
        maxN = nSamples;
    else
        maxN = min(maxN, nSamples);
    end

    if maxN<nSamples
        ids = randsample(nSamples, maxN);
    else
        ids = 1:nSamples;
    end
    
    % Number of inputs [(c), W and period]
    genes = obj.Gen.GetGenes(results(ids(1)).seq, obj.sample_genes);
    n_in = length(genes)+1;
    genes = obj.Gen.GetGenes(results(ids(1)).seq, obj.target_genes);
    n_out = length(genes); % Number of outputs
       
    samples = zeros(n_in, maxN);
    targets = zeros(n_out, maxN);
    for i = 1:maxN
        sample = results(ids(i));
        period = periods(ids(i));
        
        % Build sample and target vectors
        samples(:,i) = obj.getNNin(sample.seq, 1/period);
%         samples(:,i) = obj.getNNin(sample.seq, period);
        targets(:,i) = obj.Gen.GetGenes(sample.seq, obj.target_genes);
    end

% Normalize samples
normParams = zeros(size(samples, 1), 2);
for i = 1:size(samples, 1)
    feat = samples(i, :);
    normParams(i, :) = [mean(feat), std(feat)];
    samples(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
end
end

