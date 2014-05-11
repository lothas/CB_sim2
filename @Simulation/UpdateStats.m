function [ sim ] = UpdateStats( sim, T, X ) %#ok<INUSL>
% Update simulation statistics
Xmod = X(sim.ModCo);
Xcon = X(sim.ConCo); %#ok<NASGU>
[xNS,yNS]=sim.Mod.GetPos(Xmod,'NS'); %#ok<NASGU>
    
% Count step
sim.StepsTaken = sim.StepsTaken+1;

% Update max/min slope only if the robot is really
% walking not if it "leaped" while falling
if abs(Xmod(2)-Xmod(1))<0.5*pi % leg aperture < 90 degrees
    % Get slope on spot of "hind" leg
    Slope = sim.Env.SurfSlope(min(sim.Mod.xS,xNS))*180/pi;
    if Slope>sim.MaxSlope
        sim.MaxSlope = Slope;
    end
    if Slope<sim.MinSlope
        sim.MinSlope = Slope;
    end
end

% Check "end-game" conditions
if abs(Xmod(2)-Xmod(1))<0.00001
    sim.Out.Type = 2;
    sim.Out.Text = 'Step length too small';
    sim.StopSim = 1;
end

if sim.Mod.LegShift>0
    % Robot hit the ground before extending the leg
    sim.Out.Type = 3;
    sim.Out.Text = 'Robot hit the ground before extending the leg';
    sim.StopSim = 1;
end

% Check slope that the robot is walking on
slope = sim.Env.SurfSlope(sim.Mod.xS);
if abs(slope-sim.Env.end_slope)<0.0001 && isempty(sim.Steps2Slope)
    if sim.Env.start_slope == sim.Env.end_slope
        sim.Steps2Slope = sim.StepsTaken - 1;
    else
        sim.Steps2Slope = sim.StepsTaken;
    end
end

if sim.EndCond(1) == 1 && ~isempty(sim.Steps2Slope)
    if sim.StepsTaken-sim.Steps2Slope >= sim.EndCond(2)
        sim.Out.Type = 4;
        sim.Out.Text = ['Finished taking ',num2str(sim.StepsTaken),...
            ' steps on slope of ',num2str(sim.Env.end_slope),' rad'];
        sim.StopSim = 1;
    end        
end

end

