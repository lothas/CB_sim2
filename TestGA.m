function [  ] = TestGA(  )
GA = MOOGA(10,500);
% GA.FileIn = 'TestGA_05_25_16_32.mat';
GA.FileOut = ['TestGA_',datestr(now,'mm_dd_hh_MM'),'.mat'];
GA.Graphics = 0;
GenType = 2;

% Set up the genome
switch GenType
    case 1
        % Pulses controller (only 1 pulse to hip)
        Keys = {'omega0','P_LegE','Pulses';...
                       1,       1,   [1,2]};
        Range = {0.5, 0.55, [-20, 0, 0.01]; % Min
                   2, 0.85, [20, 0.99, 0.99]}; % Max
    case 2
        % Controller with push-off, swing pulse + limited ankle pulse
        Keys = {'omega0','P_LegE','ExtPulses','Pulses','Pulses';...
                       1,       1,      [1,1],   [1,1],   [1,2]};
        Range = {0.5, 0.55, [-30, 0.005], [-2, 0, 0.01], [0, 0, 0.01]; % Min
                   2, 0.85, [0, 0.005], [2, 0.99, 0.99], [20, 0.99, 0.99]}; % Max
end
GA.Gen = Genome(Keys, Range);

% Set up the simulation
GA.Sim = Simulation();
GA.Sim.Graphics = GA.Graphics;
GA.Sim.EndCond = 2; % Run until converge

% Set up the compass biped model
GA.Sim.Mod = GA.Sim.Mod.Set('damp',0,'I',0);

% Set up the terrain
start_slope = 0;
GA.Sim.Env = GA.Sim.Env.Set('Type','inc','start_slope',start_slope);

% Initialize the controller
GA.Sim.Con = GA.Sim.Con.ClearTorques();
GA.Sim.Con.FBType = 0;
GA.Sim.IC = [start_slope, start_slope, 0, 0, 0];

% Simulation parameters
GA.Sim = GA.Sim.SetTime(0,0.05,30);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;
% Sim.Con = Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));
GA.Sim.Con = GA.Sim.Con.HandleExtFB(GA.Sim.IC(GA.Sim.ModCo),...
                                    GA.Sim.IC(GA.Sim.ConCo));
                                
% Fitness functions
GA.NFit = 3;
GA.FitFcn = {@GA.HeightFit;
             @GA.VelFit;
             @GA.NrgEffFit};

GA = GA.InitGen();
GA = GA.Run();  
end

