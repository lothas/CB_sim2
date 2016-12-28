function [id_corr, id_conv, id_per, periods] = ...
    getConverged(obj, results, nSims)
%GETCONVERGED Gets the IDs of simulations that converged in general and
%within the desired period range in particular

    perError1_1_thresh = 0.1;
    perError1_2_thresh = 1e-6;
    perError2_1_thresh = 0.0001;

    detOsc = false(nSims,1); % Oscillations detected
    perOK1 = false(nSims,1); % Period calculation error 1
    perOK2 = false(nSims,1); % Period calculation error 2
    perOK3 = false(nSims,1); % Periodic neurons
    for i = 1:nSims
        detOsc(i) = all(~isnan(results(i).periods));
        if detOsc(i)
            perOK1(i) = max(results(i).perError1)<perError1_1_thresh;
            perOK2(i) = max(results(i).perError2)<perError2_1_thresh;
        else
            % No period was detected
            perOK1(i) = max(results(i).perError1)<perError1_2_thresh;
            perOK2(i) = all(results(i).perOK2);
        end

        oscN = results(i).neuronOsc;
        perOK3(i) = true;
        % Good combinations
        for j = 1:obj.nNeurons/2
            % If these pairs are non-oscillatory, their corresponding output
            % will be stationary
            perOK3(i) = perOK3(i) & (oscN(2*(j-1)+1) | oscN(2*j));
        end
    end
    
    % We define good results as those confirmed by both methods,
    % correctly detected (period detected and neurons are oscillatory),
    % and all signals have the same period
    periods = horzcat(results.periods);
    id_corr = find(perOK1 & perOK2 & ...
        (detOsc & perOK3 | ~detOsc & ~perOK3) & ...
        (~detOsc | abs(diff(periods))' < 0.5*min(periods(:))));
    results = results(id_corr);

    periods = horzcat(results.periods);
    % speriods = sortrows(periods,-1);
    % plot(speriods(:,1),speriods(:,2))
    periods = mean(periods)';
    id_conv = find(detOsc(id_corr));
    id_per = find(periods >= obj.perLimOut(1) & periods <= obj.perLimOut(2));

%     % How many converged to a signal of a certain period?
%     all_periods = [results.periods];
%     periods = max(all_periods);
%     % Remove ones that didn't converge
%     periods(any(isnan(all_periods))) = NaN;
%     id_conv = find(~isnan(periods));
%     nconv = length(id_conv);
%     
%     % Remove ones where we didn't detect a unique period
%     mPer = min(all_periods);
%     MPer = max(all_periods);
%     periods_diff = min(mod(MPer-mPer, mPer), mod(mPer-MPer, mPer));
%     periods(periods_diff > 0.1*periods) = NaN;
%     % Remove ones with too high/low period
%     periods(periods<0.05*obj.perLimOut(1) | ...
%         periods>20*obj.perLimOut(2)) = NaN;
% 
%     % Get CPGs with periods within range
%     id_per = find(periods>=obj.perLimOut(1) ...
%         & periods<=obj.perLimOut(2));    

    disp(['Correctly detected ',int2str(length(results)), ' out of ', ...
        int2str(nSims), ' simulations (', ...
        num2str(length(results)/nSims*100,'%.1f'),'%)']);
    disp(['Converged ',num2str(length(id_conv)), ' out of ', ...
        int2str(length(results)), ' (', ...
        num2str(length(id_conv)/length(results)*100), '%)']);
    disp(['CPGs in period range: ',num2str(length(id_per)), ' out of ', ...
        int2str(length(results)), ' (', ...
        num2str(length(id_per)/length(results)*100), '%)'])
end

