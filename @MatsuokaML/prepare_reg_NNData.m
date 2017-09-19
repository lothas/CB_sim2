function [ samples, targets, normParams ] = ...
    prepare_reg_NNData(obj, CPG_Type, filenames, maxN)
%PREPARENNDATA Prepares data to train a neural network using stored results

% *) 'CPG_Type' - can be '2N_CPG' or '4N_CPG'
%           in '2N_CPG' case the ankle period will always be NaN and we
%           dont need to care about it.

    % Load data from files
    data = load(filenames{1},'results');
    
    % Keep only results that converged
    switch CPG_Type
        case '2N_CPG'
            new_periods = horzcat(data.results(:).periods);
            new_periods = new_periods(2,:); % take only hip periods
            osc_ids = ~isnan(new_periods);

            good_ids = osc_ids;
        case '4N_CPG'
        %     new_periods = data.periods;
            new_periods = horzcat(data.results(:).periods);

            % Filter CPG's where not both signals oscillating:
            osc_ids = ~isnan(new_periods);
            osc_ids = osc_ids(1,:) & osc_ids(2,:);

            % Filter CPG's where the is a big difference between hip and ankle:
            periods_ratios = (new_periods(1,:)./new_periods(2,:));
            diff_ids = (periods_ratios >  0.85) & (periods_ratios <  1.15); 

            good_ids = osc_ids & diff_ids;
    end
    
    results = data.results(good_ids);
    periods = new_periods(:,good_ids);
    periods = mean(periods,1);
     
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
%         samples(:,i) = obj.getNNin(sample.seq, 1/period);
        samples(:,i) = obj.getNNin(sample.seq, period);
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

