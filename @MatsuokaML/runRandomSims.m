function runRandomSims(obj, nSims, filename)
%RUNRANDOMSIMS Runs simulations of the CPG using random parameters
if exist(filename, 'file') ~= 2
    simFuncHandle = @obj.runRandomSim;
    [results(nSims), ~] = simFuncHandle();

    tic
    parfor i = 1:nSims-1
        [results(i), ~] = simFuncHandle();
    end
    toc
    
    [id_conv, id_per] = getConverged(obj, results, nSims); %#ok<NASGU,ASGLU>
    
    save(filename,'results', 'nSims', 'id_conv', 'id_per');
else
    disp(['Random results already exist in ',filename])
end

end

