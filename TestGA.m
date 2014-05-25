function [  ] = TestGA(  )
GA = MOOGA(5,50);
GA.FileOut = 'TestGA.mat';

% Set up the genome
% Pulses controller (only 1 pulse to hip)
Keys = {'omega0','P_LegE','Pulses';...
               1,       1,   [1,2]};
Range = {0.5, 0.55, [-20, 0, 0.01]; % Min
           2, 0.85, [20, 0.99, 0.99]}; % Max
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
GA.NFit = 2;
GA.FitFcn = {@GA.HeightFit;
             @GA.VelFit};
%              @GA.NrgEffFit};

GA = GA.InitGen();
GA = GA.Run();  
end

