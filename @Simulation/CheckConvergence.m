function [ sim ] = CheckConvergence( sim )
% This function looks at the initial conditions stored
% and compares them to check if the system converged
% to a limit cycle

    nPer = length(sim.stepsSS); % number of periods to check
    Converge = zeros(1,nPer); % = 1 when IC are converging
    Checked = zeros(1,nPer); % = 1 when period p has been checked
    for p = 1:nPer
        if Checked(p)
            continue;
        end
        diff = norm(sim.ICstore(:,1) - sim.ICstore(:,1+p));
        
        if diff < sim.minDiff
            % Set every period multiple of p as converging
            % e.g. for 1: 1, 2, 3, ...
            %      for 2: 2, 4, 6, ...
            Converge(p:p:end) = Converge(p:p:end) + 1;
            Checked(p:p:end) = 1;
        else
            Converge(p) = 0;
            Checked(p) = 1;
        end
    end
        
    % Slow down sim for convergence if tstep_small is set
    if sum(Converge)>0 && ~isempty(sim.tstep_small)
        sim.tstep = sim.tstep_small;
    end
    
    sim.stepsSS = sim.stepsSS + Converge;
    
    % If the sim is set to stop after convergence:
    if sim.EndCond == 2        
        Period = find(sim.stepsSS >= sim.stepsReq);
        Increasing = find(Converge == 1, 1, 'first');
        if ~isempty(Period)
            if Period(1)<=Increasing
                % shortest period reached the required num. of steps
                sim.Out.Type = 5;
                sim.Out.Text = ['Reached steady state limit cycle of period ', ...
                    num2str(Period(1)),' after ',num2str(sim.StepsTaken),' steps'];
                
                % Prepare data for Poincare computation
                sim.IClimCyc = sim.ICstore(:,1);
                sim.Period = Period(1);
                
                sim.StopSim = 1;
            % else, if a lower period is still converging keep going
            end
        end
    end
end

