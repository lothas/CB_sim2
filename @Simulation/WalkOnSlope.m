function sim = WalkOnSlope(sim, start_s, end_s, leadway, max_t)
    % WALKONSLOPE Runs a simulation of the robot walking on a slope
    %   This method runs a simulation starting from start_s slope, onto
    %   end_s slope and keeps running until convergence.
    %   The simulation takes its initial conditions from sim.IClimCyc.
    
    % Simulation parameters
    sim = sim.SetTime(0,0.15,max_t);
    sim.EndCond = 2; % Run until converge

    % Some more simulation initialization
    sim.Mod.LegShift = sim.Mod.Clearance;
    %sim.Mod.xS = 0; sim.Mod.yS = 0;
    sim.Con.FBType = 2;
    sim = sim.Init();
    sim.Env = sim.Env.Set('Type','finite','start_slope',start_s,...
                          'start_x',sim.Mod.xS+leadway,...
                          'start_y',sim.Mod.yS+tand(start_s)*leadway,...
                          'end_slope',end_s,'parK',0.03);

    sim.IC = sim.IClimCyc;

    sim.Con = sim.Con.Reset(sim.IC(sim.ConCo));
    sim.Con.lastPhi = start_s*pi/180;
    if all(sim.IC) == 0
        sim.Con = sim.Con.HandleEvent(1, sim.IC(sim.ConCo));
    else
        sim.Con = sim.Con.HandleExtFB(sim.IC(sim.ModCo),...
            sim.IC(sim.ConCo),sim.Env.SurfSlope(sim.Mod.xS));
    end

    % Simulate
    sim = sim.Run();
end

