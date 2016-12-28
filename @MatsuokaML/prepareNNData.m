function [ samples, targets, normParams ] = ...
    prepareNNData(obj, input, maxN)
%PREPARENNDATA Prepares data to train a neural network using stored results
    % Load data from files
    if ischar(input{1})
        data = load(input{1});
    else
        data = input{1};
    end
    
    new_results = data.results(data.id_corr);
    new_periods = data.periods;
    % Keep only results that converged
    results = new_results(~isnan(new_periods));
    periods = new_periods(~isnan(new_periods));
    for i = 2:numel(input)
        if ischar(input{i})
            data = load(input{i});
        else
            data = input{i};
        end
    
        new_results = data.results(data.id_corr);
        new_periods = data.periods;
        % Keep only results that converged
        results = [results, new_results(~isnan(new_periods))]; %#ok<AGROW>
        periods = [periods; new_periods(~isnan(new_periods))]; %#ok<AGROW>
    end
    nSamples = numel(results);
    
    if ismember('weights',obj.sample_genes)
        pSeq = obj.permSeq(results(1).seq);
        permMult = size(pSeq,1);
    else
        permMult = 1;
    end
    
    if nargin < 3
        maxN = permMult*nSamples;
    else
        maxN = min(maxN, permMult*nSamples);
    end

    if maxN<permMult*nSamples
        ids = randsample(nSamples, ceil(maxN/permMult));
    else
        ids = 1:nSamples;
    end
    
    % Number of inputs [(c), W and period]
    genes = obj.Gen.GetGenes(results(ids(1)).seq, obj.sample_genes);
    n_in = length(genes)+1;
    genes = obj.Gen.GetGenes(results(ids(1)).seq, obj.target_genes);
    n_out = length(genes); % Number of outputs
       
    obj.normParams = [];
    
    samples = zeros(n_in, maxN);
    targets = zeros(n_out, maxN);
    for i = 1:length(ids)
        sample = results(ids(i));
        period = periods(ids(i));
        
        % Build sample and target vectors
        pSeq = obj.permSeq(sample.seq);
        for j = 1:permMult
            samples(:,(i-1)*permMult+j) = ...
                obj.getNNin(pSeq(j,:), 1/period);
    %         samples(:,i) = obj.getNNin(sample.seq, period);
            targets(:,(i-1)*permMult+j) = ...
                obj.Gen.GetGenes(pSeq(j,:), obj.target_genes);
        end
    end

% Normalize samples
normParams = zeros(size(samples, 1), 2);
for i = 1:size(samples, 1)
    feat = samples(i, :);
    normParams(i, :) = [mean(feat), std(feat)];
    samples(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
end
end

