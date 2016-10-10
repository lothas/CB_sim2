function [results, signal] = runScaledSim(obj, inputData, inputPeriod)
%RUNSCALEDSIM Runs a simulation of the CPG using new tau parameters
% rescaled to obtain a desired period

    % Setup tau_r, tau_a, c, W and feedback gains using genome
    % and beta
    seq = inputData.seq;
%     beta = inputData.b;
    
    % Select new random period within desired range
    des_period = obj.perLim(1) + rand()*(obj.perLim(2)-obj.perLim(1));

    % Scale Tr, Ta to obtain desired period
    ratio = des_period/inputPeriod;
    seq(1) = seq(1)*ratio;
    if seq(1) < obj.Gen.Range(1,1) || seq(1) > obj.Gen.Range(2,1)
        warning('Genetic sequence out of bounds, using bounded tau gene')
        % Bound tau gene
        seq(1) = min(max(seq(1), obj.Gen.Range(1,1)), obj.Gen.Range(2,1));
    end
    
    % Run simulation
    [out, sim, signal] = obj.runSim(seq);
%     [out, sim, signal] = obj.runSim(seq, beta);

    % Prepare output:
    % Parameters
    results.seq = seq;
%     results.b = beta;
    results.b = sim.Con.beta;
    results.c = sim.Con.Amp0;
    results.Worig = sim.Con.wex;
    results.W = sim.Con.W;
    results.Tr = sim.Con.tau;
    results.Ta = sim.Con.tav;
    results.x0 = out.x0;

    % Results
    results.periods = out.periods;
    results.period_Rea = out.period_Rea;
    results.amp = out.amp;
    results.pos_work = out.pos_work;
    results.neg_work = out.neg_work;
end