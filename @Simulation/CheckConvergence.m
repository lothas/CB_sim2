function [ sim ] = CheckConvergence( sim )
% This function looks at the initial conditions stored
% and compares them to check if the system converged
% to a limit cycle

    nPer = length(sim.stepsSS); % number of periods to check
    Converge = zeros(1,nPer); % = 1 when IC are converging
    Checked = zeros(1,nPer); % = 1 when period p has been checked
    mMChanged = [0,0]; % min/Max diff was changed
    thisminMax = [1,0];
    for p = 1:nPer
        if Checked(p)
            continue;
        end
        vdiff = [sim.ICstore(sim.ModCo,1) - sim.ICstore(sim.ModCo,1+p);
                 sim.Con.PhaseDiff(sim.ICstore(sim.ConCo,1),...
                                   sim.ICstore(sim.ConCo,1+p))];
        diff = norm(vdiff);
        sim.ICdiff(p) = diff;
        
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
        
        % Check convergence progression
        if sim.doGoNoGo>0
            if diff < sim.minMaxDiff(1)
                sim.minMaxDiff(1) = diff;
                mMChanged(1) = 1;
            end
            if diff > sim.minMaxDiff(2)
                sim.minMaxDiff(2) = diff;
                mMChanged(2) = 1;
            end
            % Look for the min/Max diff in this iteration
            if diff < thisminMax(1)
                thisminMax(1) = diff;
            end
            if diff > thisminMax(2)
                thisminMax(2) = diff;
            end
        end                
    end
    
    % Check convergence progression
    if sim.doGoNoGo>0
%         if ~xor(mMChanged(1),mMChanged(2))
%             % Situation unclear
%             sim.ConvProgr = [0,0];
%             sim.minMaxDiff = thisminMax;
%         else    
        if xor(mMChanged(1),mMChanged(2))
            if mMChanged(1)
                % Simulation is converging
                sim.ConvProgr(1) = sim.ConvProgr(1) + 1;
                sim.ConvProgr(2) = 0;
                sim.minMaxDiff(2) = thisminMax(2);
                
                if sim.ConvProgr(1) >= sim.GNGThresh(1)
                    if sim.doGoNoGo == 1 && sim.infTime == 0
                        % Extend the simulation by 2 step periods
                        sim.tend = min(sim.tend + 2/sim.Con.omega,...
                                       3*sim.Out.Tend); % Don't extend forever
    %                     disp(['TIME EXTENSION! ',num2str(sim.tend,'%.1f')]);
                    end
                    if sim.doGoNoGo == 2
                        % Finish the simulation
                        sim.Out.Type = 6;
                        sim.Out.Text = 'GO: Sim was converging';
                        sim.IClimCyc = sim.ICstore(:,1);
                        sim.StopSim = 1;
                    end
                end
            end
            if mMChanged(2)
                % Simulation is diverging
                sim.ConvProgr(2) = sim.ConvProgr(2) + 1;
                sim.ConvProgr(1) = 0;
                sim.minMaxDiff(1) = thisminMax(1);
                
                if sim.ConvProgr(2) >= sim.GNGThresh(2)
                    if sim.doGoNoGo == 2
                        % Stop the simulation
                        sim.Out.Type = 7;
                        sim.Out.Text = 'NO-GO: Simulation was diverging';
%                         disp(sim.Out.Text);
                        sim.tend = sim.Out.T(end);
                    end
                end
            end
        end
        
    end

    % Slow down sim for convergence if tstep_small is set
    if sum(Converge)>0 && ~isempty(sim.tstep_small)
        sim.tstep = sim.tstep_small;
    end
    
    sim.stepsSS = sim.stepsSS + Converge;
    
    % If the sim is set to stop after convergence:
    if sim.EndCond == 2        
        Period = find(sim.stepsSS >= sim.stepsReq, 1, 'first');
        Increasing = find(Converge == 1, 1, 'first');
        if ~isempty(Period)
            if Period<=Increasing
                % shortest period reached the required num. of steps
                sim.Out.Type = 5;
                sim.Out.Text = ['Reached steady state limit cycle of period ', ...
                    num2str(Period),' after ',num2str(sim.StepsTaken),' steps'];
                
                % Prepare data for Poincare computation
                sim.IClimCyc = sim.ICstore(:,1);
                sim.Period = Period;
                
                sim.StopSim = 1;
            % else, if a lower period is still converging keep going
            end
        end
    end
end

