function sim = WalkAtSpeed(sim, s_in, max_t)
    % WALKATSPPED Runs a simulation of the robot walking with a high-level
    % set speed
    %   This method runs a simulation with the given high-level walking
    %   speed command 's_in'.
    %   The simulation takes its initial conditions from sim.IClimCyc.
    
    % Simulation parameters
    sim = sim.SetTime(0,0.15,max_t);
    sim.EndCond = 2; % Run until converge

    % Some more simulation initialization
    sim.Mod.LegShift = sim.Mod.Clearance;
    sim = sim.Init();

    sim.IC = sim.IClimCyc;

    % Apply higher-level velocity signal
    sim.Con.s_in = s_in;
    sim.Con = sim.Con.Adaptation(0);

    sim.Con = sim.Con.Reset(sim.IC(sim.ConCo));
    if all(sim.IC) == 0
        sim.Con = sim.Con.HandleEvent(1, sim.IC(sim.ConCo));
    else
        sim.Con = sim.Con.HandleExtFB(sim.IC(sim.ModCo),...
            sim.IC(sim.ConCo),sim.Env.SurfSlope(sim.Mod.xS));
    end

    % Simulate
    sim = sim.Run();
end

