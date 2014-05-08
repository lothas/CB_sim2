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
    sim.End.Type = 2;
    sim.End.Text = 'Step length too small';
    sim.StopSim = 1;
end

if sim.Mod.LegShift>0
    % Robot hit the ground before extending the leg
    sim.End.Type = 3;
    sim.End.Text = 'Robot hit the ground before extending the leg';
    sim.StopSim = 1;
end

end

