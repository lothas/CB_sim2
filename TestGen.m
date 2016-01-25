function [Sim] = TestGen()
% close all

GenType = 8;

Sim = Simulation();
Sim.Graphics = 1;
Sim.EndCond = 2; % Run until converge
tend = 60;
Sim = Sim.SetTime(0,0.05,tend);

% Set up the compass biped model
Sim.Mod = Sim.Mod.Set('damp',0,'I',0);

% Set up the terrain
start_slope = 0;
Sim.Env = Sim.Env.Set('Type','inc','start_slope',start_slope);

% Set up the controller using a genome
if GenType<8
    Sim.Con = Sim.Con.ClearTorques();
end
Sim.Con.FBType = 0;

Sim.IC = [start_slope, start_slope, 0, 0, 0];

KeyLength = [];

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
    case 4 % Impulsive controller CL
        Keys = {'IC','omega0','P_LegE','AngVelImp';...
                   4,       1,       1,          2};
        Range = {[0,-4,-4,0], 0.5, 0.55, [-4 -4]; % Min
                 [1.79,4,4,0.99], 2, 0.85, [4, 4]}; % Max
               
        alpha = 0.100952073/10; phi_0 = 0.7759402;
        theta_dot = [ -0.4640/10, -0.5330/10];
        
        %                  Init Cond         omega  P_LegE  AngVelImp
        Sequence = [alpha, theta_dot, phi_0, 1.3333,  0.61,   theta_dot];
        Sim.Con.FBImpulse = 2;
        Sim.Con.ExtP_reset = phi_0;
    case 5 % Genetic Algorithm
        Keys = {'omega0','P_LegE','ExtPulses','Pulses','Pulses';
                       1,       1,      [1,1],   [1,1],   [1,2]};
        Range = {0.5, 0.55, [-30, 0.005], [-2, 0, 0.01], [-20, 0, 0.01]; % Min
                   2, 0.85, [30, 0.005], [2, 0.99, 0.99], [20, 0.99, 0.99]}; % Max
        Sequence = [1.38444, 0.787437, 3.05554, 0.005, -1.65765, 0.463841, 0.0461248,...
                    15.4709, 0.0228823, 0.0860777];
        Sim.PMFull = 1;
    case 6 % FB Genetic Algorithm
        Keys = {'omega0','P_LegE','ExtPulses','Pulses','Pulses',...
                'kOmega_u','kOmega_d','kTorques_u','kTorques_d';...
                       1,       1,      [1,1],   [1,1],   [1,2],...
                         1,         1,           1,           1};
        Range = {0.5, 0.55, [-400, 0.005], [-100, 0, 0.01], [-400, 0, 0.01],...
                  -200, -200, [-800,-400,-800], [-800,-400,-800]; % Min
                   2, 0.85, [400, 0.005], [100, 0.99, 0.99], [400, 0.99, 0.99],...
                   200,  200, [800,400,800], [800,400,800]}; % Max
        Sequence = [1.36819, 0.632062, 0.0935159, 0.005, -2.61728,...
                    0.198314, 0.371793, 66.9988, 0.0593925, 0.0207505,...
                    1.27738, -46.3519, 105.392, -370.417, 543.179,...
                    324.62, -366.454, -33.2198];
        Sequence = [1.4213, 0.661001, 54.3419, 0.005, 0.785008, 0.153893, 0.389272,...
63.3079, 0.00520723, 0.0300207, 1.21386, 3.40628, -122.345, -387.116, 414.448,...
182.592, -375.517, 157.759];
        Sim.PMFull = 1;
        Sim.Env = Sim.Env.Set('Type','finite','start_slope',start_slope,...
                                'end_slope',15);
        Gen = Genome(Keys, Range);
        KeyLength = Gen.KeyLength;
        KeyLength.kTorques_u = 3;
        KeyLength.kTorques_d = 3;
        Sim.IC = [0.1775,  -0.1775,  -0.7517,  -0.4683,  0.9804];
        Sim.IC = [0.1300 -0.1300 -0.5891 -0.3806 0.9857];
        Sim.Con.FBType = 2;
        tend = 'inf';
    case 7 % FB Genetic Algorithm (no ankle torque)
        Keys = {'omega0','P_LegE','ExtPulses','Pulses',...
                'kOmega_u','kOmega_d','kTorques_u','kTorques_d';...
                       1,       1,      [1,1],   [2,2],...
                         1,         1,           1,           1};
        Range = {0.5, 0.55, [-800, 0.005], [-400, 0, 0.01],...
                  -200, -200, [-800,-800,-800], [-800,-800,-800]; % Min
                   2, 0.85, [800, 0.005], [400, 0.99, 0.99],...
                   200,  200, [800,800,800], [800,800,800]}; % Max
        Sequence = [1.55834, 0.5935, -165.097, 0.005, -3.99613, 0.526196, 0.237336,...
26.0039, 0.0101253, 0.0589707, -34.7468, -6.625, 279.026, 185.992, 313.939,...
513.296, 571.072, -790.815];
        Sim.PMFull = 1;
        Sim.Env = Sim.Env.Set('Type','finite','start_slope',start_slope,...
                                'end_slope',-25);
        Gen = Genome(Keys, Range);
        KeyLength = Gen.KeyLength;
        KeyLength.kTorques_u = 3;
        KeyLength.kTorques_d = 3;
%         Sim.IC = [0.1775,  -0.1775,  -0.7517,  -0.4683,  0.9804];
%         Sim.IC = [0.1300 -0.1300 -0.5891 -0.3806 0.9857];
        Sim.Con.FBType = 2;
        tend = 10;
        Sim.tstep_small = 0.001;
    case 8 % Matsuoka controller
        Sim.Con = Matsuoka;
        
        nAnkle = 1; % Number of ankle torques
        nHip = 1;   % Number of hip torques
        maxAnkle = 8;   % Max ankle torque
        maxHip = 20;    % Max hip torque
        Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        N = nAnkle+nHip;
        Mwex = -10*ones(1,4*N);
        mwex = 0*Mwex;
        
        Sim.Con.nPulses = N;
        Sim.Con.stDim = 4*N;
        Sim.Con = Sim.Con.SetOutMatrix([nAnkle,nHip]);
        Sim.Con.MinSat = [-maxAnkle,-maxHip];
        Sim.Con.MaxSat = [ maxAnkle, maxHip];
        
        Keys = {'tau','tav','beta','amp','win','wex','ks_tau',  'ks_out','IC_matsuoka';
                   1 ,   1 ,    1 , 2*N ,   1 , 4*N ,      1 ,      2*N ,           0 };
        Range = {0.1 , 0.1 ,    1 , mamp, -20 , Mwex,   -1e2 , -0.1*Mamp; % Min
                   2 ,   2 ,   30 , Mamp,  -1 , mwex,    1e2 ,  0.1*Mamp}; % Max
                   
        %            tau,   tav, beta,                   amp,    win,
        Sequence = [0.3,  0.6,   15, [1.00,4.00,10.00,1.00],     -1,...
                    [-1, 0, 0,-3,-1, 0, 0,-3], -0.0001, [-0.02 0.02 0.4 0.1]];
        %                                 wex,   s_tau,               s_out
end

if ~isempty(KeyLength)
    Gen = Genome(Keys, KeyLength, Range);
else
    Gen = Genome(Keys, Range);
end
Sim = Gen.Decode(Sim, Sequence);

% Set internal parameters (state dimensions, events, etc)
Sim = Sim.Init();

% Some more simulation initialization
Sim.Mod.LegShift = Sim.Mod.Clearance;
% Sim.Con = Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));
Sim.Con = Sim.Con.HandleExtFB(Sim.IC(Sim.ModCo), ...
                              Sim.IC(Sim.ConCo), start_slope);

                          
                          
                          
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
    close all
    % Initial conditions
    IC0 = Sim.IC(Sim.ConCo);

    % tspan
    t_start = 0;
    t_end = 3.5;
    t_step = 0.03;
    t_span = t_start:t_step:t_end;

    % Run simulation
    options = odeset('MaxStep',t_step/10,'RelTol',.5e-7,'AbsTol',.5e-8);
    [TTemp,XTemp] = ode45(@Derivative,t_span,IC0,options);
    % options = odeset('MaxStep',t_step/10,'RelTol',.5e-7,'AbsTol',.5e-8,...
    %             'Events', @MO.Events);
    % [TTemp,XTemp,TE,YE,IE] = ...
    %     ode45(@MO.Derivative,t_span,IC0,options); %#ok<ASGLU>

    y = max(XTemp(:,1:2:end),0);
    figure
    N = Sim.Con.nPulses;
    for i = 1:4
        subplot(4,1,i)
        plot(TTemp, XTemp(:,N*(i-1)+1:N*i));
    end
    figure
    plot(TTemp, Sim.Con.Output(0, XTemp', 0));
    grid on

    return
end
        
        function xdot = Derivative(t,x)
            xdot = Sim.Con.Derivative(t,x);
        end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%                     
                          
                          
                          
                          
% Simulate
Sim = Sim.Run();

[fit,out] = MOOGA.VelRangeFit(Sim);
MOOGA.NrgEffFit(Sim);
MOOGA.UphillFitRun(Sim);

% Calculate eigenvalues
if Sim.Out.Type == 5
    [EigVal,EigVec] = Sim.Poincare();
    % Do some plots
    disp(EigVal);
else
    EigVal = 2*ones(4,1);
    disp(Sim.Out.Text);
end

% Display steady state initial conditions
disp(Sim.ICstore(:,1));
end