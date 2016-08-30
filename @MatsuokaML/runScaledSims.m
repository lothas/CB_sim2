function runScaledSims(obj, inputData, filename)
%RUNSCALEDSIMS Re-runs simulations that converged outside the desired
% period range, this time with scaled temporal parameters
if exist(filename, 'file') ~= 2
    nSims = numel(inputData);
    results(nSims) = inputData(1); % Initialize storage

    tic
    per_min = obj.perLim(1);
    per_rng = obj.perLim(2) - obj.perLim(1);
    simFuncHandle = @obj.sim;
    parfor i = 1:nSims
        % Select new random period within desired range
        des_period = per_min + rand()*per_rng;

        % Scale Tr, Ta to obtain desired period
        ratio = des_period/max(inputData(i).periods);
        Tr = inputData(i).Tr*ratio;
        Ta = 5*Tr;

        % Re-run simulation
        [out, ~] = simFuncHandle(inputData(i).b, inputData(i).c, ...
                                      inputData(i).W, Tr, Ta);

        % Prepare output:
        % Parameters
        results(i).a = inputData(i).a;
        results(i).b = inputData(i).b;
        results(i).c = inputData(i).c;
        results(i).Worig = inputData(i).Worig;
        results(i).W = inputData(i).W;
        results(i).Tr = Tr;
        results(i).Ta = Ta;
        results(i).x0 = out.x0;

        % Results
        results(i).periods = out.periods;
        results(i).pos_work = out.pos_work;
        results(i).neg_work = out.neg_work;
    end
    toc
    
    [id_conv, id_per] = getConverged(obj, results, nSims); %#ok<NASGU,ASGLU>
    
    save(filename, 'results', 'nSims', 'id_conv', 'id_per');
else
    disp(['Random results already exist in ',filename])
end

end

