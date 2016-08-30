function [ samples, targets, normParams ] = ...
    prepareNNData(obj, filenames, maxN)
%PREPARENNDATA Prepares data to train a neural network using stored results
    % Number of inputs (c, W and period)
    n_in = obj.nNeurons + (obj.nNeurons-1)*obj.nNeurons + 1;
    n_out = 1; % Number of outputs (Tau_r)
    
    % Load data from files, keep only results that converged
    data = load(filenames{1});
    results = data.results(data.id_conv);
    for i = 2:numel(filenames)
        data = load(filenames{i});
        results = [results, data.results(data.id_conv)]; %#ok<AGROW>
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
    
    samples = zeros(n_in, maxN);
    targets = zeros(n_out, maxN);
    for i = 1:maxN
        sample = results(ids(i));
        
        period = max(sample.periods);
        W = sample.W;
        weights = W(~logical(eye(size(W))));

        % Build sample and target vectors
        %               Tonic input  Weights   Period
        samples(:,i) = [sample.c;    weights;  period];
        %               Tau mem. pot.
        targets(:,i) = [sample.Tr];
    end
    
    % Normalize samples
    normParams = zeros(size(samples, 1), 2);
    for i = 1:size(samples, 1)
        feat = samples(i, :);
        normParams(i, :) = [mean(feat), std(feat)];
        samples(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
    end
end

