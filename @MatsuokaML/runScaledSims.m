function runScaledSims(obj, inputData, inputPeriods, filename)
%RUNSCALEDSIMS Re-runs simulations that converged outside the desired
% period range, this time with scaled temporal parameters
if exist(filename, 'file') ~= 2
    nSims = numel(inputData);
    results(nSims) = inputData(1); % Initialize storage

    % Reduce simulation time to ~20 cycles
    obj.tEnd = obj.perLimOut(2)*20;
    simFuncHandle = @obj.runScaledSim;
    
    tic
    parfor i = 1:nSims
        [results(i), ~] = simFuncHandle(inputData(i), inputPeriods(i));
    end
    toc
    
    [id_corr, id_conv, id_per, periods] = ...
        getConverged(obj, results, nSims); %#ok<ASGLU>
    
    save(filename, 'results', 'nSims', 'id_corr', ...
        'id_conv', 'id_per', 'periods');
else
    disp(['Random results already exist in ',filename])
end

end

