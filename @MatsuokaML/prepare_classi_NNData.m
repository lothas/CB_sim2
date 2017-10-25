function [ samples, targets] = ...
    prepare_classi_NNData(obj, CPG_Type, filenames, maxN)
%PREPARE_CLASSI_NNDATA 
% Prepares data to train a classification neural network using stored results

% *) 'CPG_Type' - can be '2N_CPG' or '4N_CPG'
%           in '2N_CPG' case the ankle period will always be NaN and we
%           dont need to care about it.
    
    results = [];
    periods = [];
    
    % Load data from files
    for i = 1:numel(filenames)
        data = load(filenames{i});
        new_periods = horzcat(data.results(:).periods);
        results = [results, data.results]; %#ok<AGROW>
        periods = [periods, new_periods]; %#ok<AGROW>
    end
    
    switch CPG_Type
        case '2N_CPG'
            periods = periods(2,:); % take only hip period
        case '4N_CPG'
            % Filter CPG's where not both signals oscillating:
            osc_ids = ~isnan(periods);
            osc_ids = osc_ids(1,:) & osc_ids(2,:);

            % Filter CPG's where the is a big difference between hip and ankle:
            periods_ratios = (periods(1,:)./periods(2,:));
            diff_ids = (periods_ratios >  0.85) & (periods_ratios <  1.15); 

            good_ids = osc_ids & diff_ids;
            
            periods = mean(periods,1);
            
            % classifi all of the "bad" CPGs as "n-osc":
            periods(1,~good_ids) = NaN(1,sum(~good_ids));
    end
    
    % get the number of samples:
    nSamples = numel(results);
    
    % get random 'maxN' samples from all of the samples:
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
    
    % %  Number of inputs [(c), W and period]
    genes = obj.Gen.GetGenes(results(ids(1)).seq, obj.sample_genes);
    n_in = length(genes);
        % NOTE: 'period' is NOT an input for the classification NN
    
    % there is always one NN_out which is the probability of the sample
    %   to be in each class ('n_osc' and 'osc' classes)
    n_out = 2; % Number of outputs
       
    samples = zeros(n_in, maxN);
    targets = zeros(n_out, maxN);
    
    for i = 1:maxN
        sample = results(ids(i));
        period = periods(ids(i));
        
        % Build sample and target vectors
        temp = obj.getNNin(sample.seq, period);
            % remove the periods from the inputs:
        samples(:,i) = temp(1:end-1,:);
        
        % devide to classes:
        %   targ=[1;0] if CPG is n_osc
        %   targ=[0;1] if CPG is osc 
        targets(:,i) = [isnan(period);...
            ~isnan(period)];
    end
    
%     % Normalize samples
%     normParams = zeros(size(samples, 1), 2);
%     for i = 1:size(samples, 1)
%         feat = samples(i, :);
%         normParams(i, :) = [mean(feat), std(feat)];
%         samples(i, :) = (feat - normParams(i, 1))/normParams(i, 2);
%     end

end

