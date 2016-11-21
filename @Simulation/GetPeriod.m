function Period = GetPeriod(sim, tspan)
%GETPERIOD Returns the walking period (calculates it if necessary)
if size(sim.Period,2) == 2
    Period = sim.Period;
    return
end

if nargin<2
    tspan = 1;
end

if length(tspan) == 1
    ind_end = length(sim.Out.T);
    if tspan<=1
        % Get the last "tspan" percentage of the simulation
        tspan = floor((1-tspan)*ind_end) + 1 : ind_end;
    else
        tspan = ind_end - tspan + 1 : ind_end;
    end
end
    
% Calculate walking period
T = sim.Out.T(tspan);
SuppPos = sim.Out.SuppPos(tspan,:);

StepTimes = T(diff(SuppPos(:,1))~=0);
Periods = diff(StepTimes);

% There seems to be a problem calculating
% the step period sometimes so we'll try
% another "back-up" way
if isempty(Periods)
    % run the simulation for one step
    psim = copy(sim);
    if isempty(sim.Period)
        psim.EndCond = [1,10];
    else
        psim.EndCond = [1,sim.Period(1)];
    end
    psim.IC = sim.IClimCyc;
    psim = psim.Init();
    psim = psim.Run();
    Periods = psim.Out.T(end) - psim.Out.T(1);
else    
    mean_P = mean(Periods);
    std_P = std(Periods);
    
    Periods(Periods > mean_P+std_P | Periods < mean_P-std_P) = [];
end


if size(sim.Period,2) == 1
    Period = [sim.Period, mean(Periods)];
else
    % Return period 1 as default
    Period = [1, mean(Periods)];
end


end

