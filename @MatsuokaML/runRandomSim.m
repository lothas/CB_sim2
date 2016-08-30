function [results, signal] = runRandomSim(obj)
%RUNRANDOMSIM Runs a simulation of the CPG using random parameters
    % Setup oscillators' parameters
    [a, b, c, Worig, W, Tr, Ta] = obj.getRandPar();
    
    % Run simulation
    [out, signal] = obj.sim(b, c, W, Tr, Ta);

    % Prepare output:
    % Parameters
    results.a = a;
    results.b = b;
    results.c = c;
    results.Worig = Worig;
    results.W = W;
    results.Tr = Tr;
    results.Ta = Ta;
    results.x0 = out.x0;

    % Results
    results.periods = out.periods;
    results.pos_work = out.pos_work;
    results.neg_work = out.neg_work;
end

