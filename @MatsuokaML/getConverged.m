function [id_conv, id_per, periods] = getConverged(obj, results, nSims)
%GETCONVERGED Gets the IDs of simulations that converged in general and
%within the desired period range in particular
    % How many converged to a signal of a certain period?
    all_periods = [results.periods];
    periods = max(all_periods);
    % Remove ones that didn't converge
    periods(any(isnan(all_periods))) = NaN;
    id_conv = find(~isnan(periods));
    nconv = length(id_conv);
    
    % Remove ones where we didn't detect a unique period
    periods_diff = mod(max(all_periods)-min(all_periods), ...
                        min(all_periods));
    periods(periods_diff > 0.1*periods) = NaN;
    % Remove ones with too high/low period
    periods(periods<0.05*obj.perLimOut(1) | ...
        periods>20*obj.perLimOut(2)) = NaN;
    
    disp(['Converged ',num2str(nconv), ' out of ',int2str(nSims), ...
        ', i.e. ',num2str(nconv/nSims*100), '%'])

    % Get CPGs with periods within range
    id_per = find(periods>=obj.perLimOut(1) ...
        & periods<=obj.perLimOut(2));

    disp(['CPGs in period range: ',num2str(length(id_per)), ' out of ', ...
        int2str(nSims), ', i.e. ', ...
        num2str(length(id_per)/nSims*100), '%'])
end

