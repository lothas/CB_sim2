function [Sim] = ImpulseSim(varargin)
Sim = Simulation();
Sim.Graphics = 1;
Sim.EndCond = 2; % Run until convergence

% Set up the compass biped model
Sim.Mod = Sim.Mod.Set('damp',0,'I',0); % no damping, point-mass

% Set up the terrain
start_slope = 0;
Sim.Env = Sim.Env.Set('Type','inc','start_slope',start_slope);

% Impulsive controller
phi_0 = 0.7759402;

Sim.Con = Sim.Con.ClearTorques(); % no torques
Sim.Con = Sim.Con.Set('omega0', 1.03333); % 1/T;T =0.8895

switch nargin
    case 0
        % Default test case
%         Sim.Con.FBImpulse = 2; % set ang. vel. to certain value
%         alpha = 0.110952073;
% %         theta_dot = [ -0.4640, -0.5330 ];
%         theta_dot = [ -0.5040, -0.2530 ];
%         impulse = theta_dot;
        
        Sim.Con.FBImpulse = 0; % set ang. vel. to certain value
        Sim.Con.FBType = 0;
        Sim.EndZMP = 0;
        
        T = 0.679875183484506; omega = 1/T;
        Sim.Con = Sim.Con.Set('omega0', omega);
        
        alpha = 0.08777523036753;
        theta_dot = [ -0.386077676960781, -0.359050627940161 ];
        impulse = theta_dot;
        
        delta = [-0.019877882616433  -0.126797754378412];
        delta_joint = [delta(1)+delta(2) delta(2)];
        duration = 0.05;
        amp = delta_joint/duration;
        Sim.Con = Sim.Con.AddPulse('joint',1,'amp',amp(1),'offset','ext','dur',duration);
        Sim.Con = Sim.Con.AddPulse('joint',2,'amp',amp(2),'offset','ext','dur',duration);
    case 2
        % Full impulse feedback case
        % Input: Initial alpha and desired ang. vel.
        Sim.Con.FBImpulse = 2; % set ang. vel. to certain value
        alpha = varargin{1};
        theta_dot = varargin{2};
        impulse = theta_dot;
    case 3
        % Open loop impulse case
        % Input: Initial alpha, ang. vel. and impulse value
        Sim.Con.FBImpulse = 1; % Add impulse value to ang. vel.
        alpha = varargin{1};
        theta_dot = varargin{2};
        impulse = varargin{3};
    case 4
        % Open loop quasi-impulse case
        % Input: Initial alpha, ang. vel., desired impulse value and
        %        dt = pulse duration
        Sim.Con.FBImpulse = 0; % No impulse added (only short pulses below)
        Sim.Con.FBType = 0; % No feedback
        Sim.EndZMP = 0;
        
        T = 1; omega = 1/T;
        Sim.Con = Sim.Con.Set('omega0', omega);
        
        alpha = varargin{1};
        theta_dot = varargin{2};
        impulse = varargin{3};
        delta = varargin{4};
        delta_phi = delta*omega;
        amp = impulse/delta_phi;
        Sim.Con = Sim.Con.AddPulse('joint',1,'amp',amp(1), ...
                                   'offset','ext','dur',delta_phi);
        Sim.Con = Sim.Con.AddPulse('joint',2,'amp',amp(2), ...
                                   'offset','ext','dur',delta_phi);
end

theta = [start_slope+alpha,start_slope-alpha];
Sim.Con.ExtP_reset = phi_0; % enable phase reset for CPG
Sim.Con.AngVelImp = impulse; % set CPG impulse value

% Simulation parameters
Sim = Sim.SetTime(0,0.05,60);
Sim.IC = [theta, theta_dot, phi_0];

% Set internal parameters (state dimensions, events, etc)
Sim = Sim.Init();

% Some more simulation initialization
Sim.Mod.LegShift = Sim.Mod.Clearance;
Sim.Con.HandleExtFB(Sim.IC(Sim.ModCo), Sim.IC(Sim.ConCo), start_slope);
end