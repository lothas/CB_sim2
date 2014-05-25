function [Sim] = TestGen()
% close all

GenType = 2;

Sim = Simulation();
Sim.Graphics = 1;
Sim.EndCond = 2; % Run until converge

% Set up the compass biped model
Sim.Mod = Sim.Mod.Set('damp',0,'I',0);

% Set up the terrain
start_slope = 0;
Sim.Env = Sim.Env.Set('Type','inc','start_slope',start_slope);

% Set up the controller using a genome
Sim.Con = Sim.Con.ClearTorques();
Sim.Con.FBType = 0;
Sim.IC = [start_slope, start_slope, 0, 0, 0];
switch GenType
    case 1 % Event triggered controller
        Keys = {'IC','omega0','P_LegE','ExtPulses','ExtPulses';...
                   4,       1,       1,      [1,1],      [1,2]};
        Range = {[0,-2,-2,0], 0.5, 0.55, [-20, 0.01], [-20, 0.01]; % Min
                 [0.79,2,2,0.99], 2, 0.85, [20, 0.99], [20, 0.99]}; % Max
               
        T = 0.779875183484506; alpha = 0.08777523036753; phi_0 = 0.7759402;
        theta_dot = [ -0.386077676960781, -0.359050627940161 ];
        Sim.Con.ExtP_reset = phi_0;
        delta = [-0.019877882616433  -0.126797754378412];
        delta_joint = [delta(1)+delta(2) delta(2)];
        duration = 0.05; amp = delta_joint/duration;
        
        %                  Init Cond        omega P_LegE  ExtPulses ankle    ExtPulses hip  
        Sequence = [alpha, theta_dot, phi_0, 1/T,  0.61,  amp(1), duration, amp(2), duration];
    case 2 % Pulses controller
        Keys = {'omega0','P_LegE','Pulses','Pulses';...
                       1,       1,   [1,1],   [1,2]};
        Range = {0.5, 0.55, [-20, 0, 0.01], [-20, 0, 0.01]; % Min
                   2, 0.85, [20, 0.99, 0.99], [20, 0.99, 0.99]}; % Max
                   
        %            omega  P_LegE         Pulses ankle             Pulses hip  
        Sequence = [1.2666, 0.5973, -7.3842, 0.1268, 0.07227, 5.1913, 0.1665, 0.0537];
        Sim.IC = [0,0,0,0,0];
    case 3 % Impulsive controller
        Keys = {'IC','omega0','P_LegE','AngVelImp';...
                   4,       1,       1,          2};
        Range = {[0,-2,-2,0], 0.5, 0.55, [-2 -2]; % Min
                 [0.79,2,2,0.99], 2, 0.85, [2, 2]}; % Max
               
        alpha = 0.100952073; phi_0 = 0.7759402;
        theta_dot = [ -0.4640, -0.5330 ]; delta = [-0.019877882616433  -0.126797754378412];
        
        %                  Init Cond         omega  P_LegE  AngVelImp
        Sequence = [alpha, theta_dot, phi_0, 1.3333,  0.61,   delta];
        Sim.Con.FBImpulse = 1;
        Sim.Con.ExtP_reset = phi_0;
end
Gen = Genome(Keys, Range);
Sim = Gen.Decode(Sim, Sequence);

% Simulation parameters
Sim = Sim.SetTime(0,0.05,60);

% Set internal parameters (state dimensions, events, etc)
Sim = Sim.Init();

% Some more simulation initialization
Sim.Mod.LegShift = Sim.Mod.Clearance;
% Sim.Con = Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));
Sim.Con = Sim.Con.HandleExtFB(Sim.IC(Sim.ModCo),Sim.IC(Sim.ConCo));

% Simulate
Sim = Sim.Run();

% Calculate eigenvalues
if Sim.Out.Type == 5
    [EigVal,EigVec] = Sim.Poincare();
    % Do some plots
    disp(EigVal);
else
    EigVal = 2*ones(4,1);
    disp(Sim.Out.Text);
end

end