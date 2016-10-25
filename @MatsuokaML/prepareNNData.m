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
    
%     % Normalize weights
%     weight_id = find(strcmp(obj.sample_genes,'weights'),1,'first');
%     norm_weights = 0;
%     if ~isempty(weight_id)
%         norm_weights = 1;
%         sample_genes1 = {'weights'};
%         sample_genes2 = obj.sample_genes;
%         sample_genes2(weight_id) = [];
%     end
    
    samples = zeros(n_in, maxN);
    targets = zeros(n_out, maxN);
    for i = 1:maxN
        sample = results(ids(i));
        period = periods(ids(i));
        
        % Get selected genes from genetic sequence
%         if norm_weights
%             genes1 = obj.Gen.GetGenes(sample.seq, sample_genes1); % weights
%             amp = obj.Gen.GetGenes(sample.seq, {'amp'}); % weights
%             for r = 1:obj.nNeurons
%                 range = (obj.nNeurons-1)*(r-1)+1:(obj.nNeurons-1)*r;
%                 this_amp = amp(r);
%                 amps = amp; amps(r) = [];
%                 genes1(range) = genes1(range)*this_amp./amps;
%             end
%             genes1 = min(max(genes1, -10),100); % Bound the weights
%             genes2 = obj.Gen.GetGenes(sample.seq, sample_genes2);
%             genes = [genes1, genes2];
%         else
            genes = obj.Gen.GetGenes(sample.seq, obj.sample_genes);
%         end
        target = obj.Gen.GetGenes(sample.seq, obj.target_genes);
        
        % Build sample and target vectors
        samples(:,i) = [genes';  period];
        %              Tau mem. pot.
        targets(:,i) = target;
    end
    
    % Normalize samples
    normParams = zeros(size(samples, 1), 2);
    for i = 1:size(samples, 1)
        feat = samples(i, :);
        normParams(i, :) = [mean(feat), std(feat)];
        samples(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
    end
end

