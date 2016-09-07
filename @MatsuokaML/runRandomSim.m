function [results, signal] = runRandomSim(obj)
%RUNRANDOMSIM Runs a simulation of the CPG using random parameters
    % Setup oscillators' parameters
%     [a, b, c, Worig, W, Tr, Ta] = obj.getRandPar();

    % Setup tau_r, tau_a, c, W and feedback gains using genome
    seq = obj.Gen.RandSeq();
    
    % Set random b
    if rand()>0.7
        beta = min(max(0.6+0.1*randn(),0.2),0.8);
    else
        beta = min(max(2.5+randn(),0.8),8);
    end
    
    % Run simulation
    [out, sim, signal] = obj.runSim(seq, beta);

    % Prepare output:
    % Parameters
    results.seq = seq;
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

