function [Sim] = TestSim()
close all

Sim = Simulation();
Sim.Graphics = 1;
Sim.EndCond = 2;

% Set up the compass biped model
Sim.Mod = Sim.Mod.Set('damp',0,'I',0);

% Set up the terrain
start_slope = 0;
Sim.Env = Sim.Env.Set('Type','inc','start_slope',start_slope);

% Set up the controller
% Sim.Con = Sim.Con.ClearTorques();
% Sim.Con = Sim.Con.Set('omega0',1.2666,'P_LegE',0.5973,'FBType',0);
% Sim.Con = Sim.Con.AddPulse('joint',1,'amp',-7.3842,'offset',0.1268,'dur',0.07227);
% Sim.Con = Sim.Con.AddPulse('joint',2,'amp',5.1913,'offset',0.1665,'dur',0.0537);

% Let's try an impulsive controller
Sim.Con = Sim.Con.ClearTorques();
Sim.Con = Sim.Con.Set('omega0', 1.3333,'P_LegE',0.61); % 1/T;T =0.8895
Sim.Con.FBImpulse = 2;
theta_dot = [ -0.4640, -0.5330 ];
alpha = 0.100952073;
thetta = [start_slope+alpha,start_slope-alpha];
phi_0 = 0.7759402;
Sim.Con.ExtP_reset = phi_0;
Sim.Con.AngVelImp = theta_dot;

% Simulation parameters
Sim = Sim.SetTime(0,0.05,60);
% Sim.IC = [0.13, -0.1, -0.4, -0.25, 0];
% Sim.IC = [0., 0., 0., 0., 0.];
% Sim.IC = [0.1393442, -0.1393442, -0.5933174, -0.4680616, 0.8759402];
Sim.IC = [thetta, theta_dot, phi_0];

% Set internal parameters (state dimensions, events, etc)
Sim = Sim.Init();

% Some more simulation initialization
Sim.Mod.LegShift = Sim.Mod.Clearance;
Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));

% Simulate
Sim = Sim.Run();

% Calculate eigenvalues
[EigVal,EigVec] = Sim.Poincare();

% Do some plots
disp(EigVal);
end