function [ sim ] = Init( sim )
% Initialize simulation properties
    % Set states
    sim.stDim = sim.Mod.stDim + sim.Con.stDim; % state dimension
    sim.ModCo = 1:sim.Mod.stDim; % Model coord. indices
    sim.ConCo = sim.Mod.stDim+1:sim.stDim; % Contr. coord. indices

    % Set events
    sim.nEvents = sim.Mod.nEvents + sim.Con.nEvents;
    sim.ModEv = 1:sim.Mod.nEvents; % Model events indices
    sim.ConEv = sim.Mod.nEvents+1:sim.nEvents; % Contr. events indices
    
    sim.StopSim = 0;
        
    % Set render params
    if sim.Graphics == 1
        if sim.Fig == 0
            sim.Once = 1;
        end
        
        % Init window size params
        scrsz = get(0, 'ScreenSize');
        if scrsz(3)>2*scrsz(4) % isunix()
            % If 2 screens are used in Linux
            scrsz(3) = scrsz(3)/2;
        end
        sim.FigWidth = scrsz(3)-250;
        sim.FigHeight = scrsz(4)-250;
        sim.AR = sim.FigWidth/sim.FigHeight;
        [sim.COMx0,sim.COMy0] = sim.Mod.GetPos(sim.IC(sim.ModCo),'COM');
        
        % Init world size params
        sim.FlMin = sim.COMx0-1.25*sim.AR*sim.Mod.L;
        sim.FlMax = sim.COMx0+1.25*sim.AR*sim.Mod.L;
        sim.HeightMin = sim.COMy0-1/sim.AR*sim.Mod.L;
        sim.HeightMax = sim.COMy0+4/sim.AR*sim.Mod.L;

        % Init torque display params
        sim.nOuts = length(sim.Con.Amp0);
        if sim.nOuts>0
            % Set number of steps so a whole cycle of the oscillator
            % will be included
            sim.nTsteps = ceil(sim.Con.GetPeriod()/sim.tstep);
            sim.Ttime = linspace(sim.FlMax*0.8,sim.FlMax*0.95,sim.nTsteps);
            sim.Thold = zeros(sim.nOuts,sim.nTsteps);
            sim.Tbase = (sim.HeightMax+sim.HeightMin)/2;
            sim.Tscale = 0.1*(sim.HeightMax-sim.HeightMin)/max(abs(sim.Con.Amp0));
            sim.hTorques = zeros(sim.nOuts,1);
        end
    end
    
    sim.Mod.curSpeed = 'Computing...';
    sim.StepsTaken = 0;
    
    % Adapt CPG (if adaptive)
    sim.Con = sim.Con.Adaptation([0 2*sim.Env.SurfSlope(sim.Mod.xS) 0 0 0]);
end

