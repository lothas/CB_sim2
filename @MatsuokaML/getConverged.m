function [id_conv, id_per] = getConverged(obj, results, nSims)
%GETCONVERGED Gets the IDs of simulations that converged in general and
%within the desired period range in particular
    % How many converged to a signal of a certain period?
    id_conv = find(~any(isnan([results.periods])));
    nconv = length(id_conv);

    disp(['Converged ',num2str(nconv), ' out of ',int2str(nSims), ...
        ', i.e. ',num2str(nconv/nSims*100), '%'])

    % Get CPGs with periods within range
    mean_periods = mean([results.periods])';
    id_per = find(mean_periods>=obj.perLimOut(1) ...
        & mean_periods<=obj.perLimOut(2));

    disp(['CPGs in period range: ',num2str(length(id_per)), ' out of ', ...
        int2str(nSims), ', i.e. ', ...
        num2str(length(id_per)/nSims*100), '%'])
end

