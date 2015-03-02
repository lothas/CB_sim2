function [  ] = TestGApassive(  )
GA = MOOGA(5,250);
GA.FileIn = 'TestGAp_03_02_13_20.mat';
GA.FileOut = ['TestGAp_',datestr(now,'mm_dd_hh_MM'),'.mat'];
GA.Graphics = 0;
GA.ReDo = 1;

% Set up the genome
Keys = {'IC','omega0';...
           4,       1};
Range = {[0,-2,-2.5,0.6], 1.2; % Min
         [0.5,0,0,0.99], 1.4}; % Max
GA.Gen = Genome(Keys, Range);
GA.JOAT = 0;

% Set up the simulation
GA.Sim = Simulation();
GA.Sim.Graphics = GA.Graphics;
GA.Sim.EndCond = 2; % Run until converge

% Set up the compass biped model
GA.Sim.Mod = GA.Sim.Mod.Set('damp',0,'I',0);

% Set up the terrain
start_slope = -1;
GA.Sim.Env = GA.Sim.Env.Set('Type','inc','start_slope',start_slope);

% Initialize the controller
GA.Sim.Con = GA.Sim.Con.ClearTorques();
GA.Sim.Con.FBType = 0;
GA.Sim.Con.ExtP_reset = 0;
GA.Sim.IC = [start_slope*pi/180, start_slope*pi/180, 0, 0, 0];

% Simulation parameters
GA.Sim = GA.Sim.SetTime(0,0.05,30);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;
% Sim.Con = Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));
% GA.Sim.Con = GA.Sim.Con.HandleExtFB(GA.Sim.IC(GA.Sim.ModCo),...
%                                     GA.Sim.IC(GA.Sim.ConCo));
                                
% Fitness functions
GA.FitFcn = {1, @MOOGA.HeightFit;
             2, @MOOGA.VelFit};
GA.FitIDs = [1,2];
GA.NFit = size(GA.FitFcn,1);

GA = GA.InitGen();
GA = GA.Run();  
end

