function GenFBtuner()
% close all
%           omega    P_LegE   pulse [Amp, Offset, Duration]
Sequence = [1.47469, 0.727897,...
            -24.2109, 0.539044, 0.201354,... % Ankle pulse 1
            8.91417, 0.532459, 0.259763,... % Ankle pulse 2
            4.85994, 0.591071, 0.15863,... % Ankle pulse 3
            40.0000, 0.058046, 0.0247,... % Hip pulse 1
                  6, 0.162751, 0.17404,... % Hip pulse 2
           -6.8021, 0.560096, 0.297654,... % Hip pulse 3
            3,3,... % omega gain up/down
            -330,-100,-100,50,50,0,... % pulse gains up
            0,0,0,0,0,0]; % pulse gains down
        
% Sequence = [1.47469, 0.727897,...
%             -24.2109, 0.539044, 0.201354,... % Ankle pulse 1
%             8.91417, 0.532459, 0.259763,... % Ankle pulse 2
%             4.85994, 0.591071, 0.15863,... % Ankle pulse 3
%             56.5475, 0.438046, 0.193643,... % Hip pulse 1
%            -50.3699, 0.462751, 0.17404,... % Hip pulse 2
%            -13.8021, 0.560096, 0.297654,... % Hip pulse 3
%             0,0,... % omega gain up/down
%             0,0,0,0,0,0,... % pulse gains up
%             0,0,0,0,0,0]; % pulse gains down

Sim = Simulation();
Sim.Graphics = 1;
Sim.EndCond = 2; % Run until converge
tend = 'inf';

% Set up the compass biped model
Sim.Mod = Sim.Mod.Set('damp',0,'I',0);

% Set up the terrain
start_slope = 0;
Sim.Env = Sim.Env.Set('Type','inc','start_slope',start_slope);

% Simulation parameters
Sim.IC = [start_slope, start_slope, 0, 0, 0];
% Sim.IC = [0.186672, -0.186672, -0.822393, -0.421087, 0.941741];

% Set up the genome
% Controller with swing pulse + limited ankle pulse and feedback
NAnkleT = 3;
NHipT = 3;
MaxAnkleT = 150;
MaxHipT = 150;
TorqueFBMin = [-1000*ones(1,NAnkleT),-1000*ones(1,NHipT)];
TorqueFBMax = -TorqueFBMin;
Keys = {'omega0','P_LegE',    'Pulses',    'Pulses',...
    'kOmega_u','kOmega_d','kTorques_u','kTorques_d';...
               1,       1, [NAnkleT,1],   [NHipT,2],...
             1,         1,           1,           1};
Range = {0.7, 0.55, [-MaxAnkleT, 0, 0.01], [-MaxHipT, 0, 0.01],...
    -3, -3, TorqueFBMin, TorqueFBMin; % Min
         1.5, 0.85, [MaxAnkleT, 0.99, 0.99], [MaxHipT, 0.99, 0.99],...
    3,  3, TorqueFBMax, TorqueFBMax}; % Max

Gen = Genome(Keys, Range);
KeyLength = Gen.KeyLength;
KeyLength.kTorques_u = length(TorqueFBMin);
KeyLength.kTorques_d = length(TorqueFBMin);
Gen = Genome(Keys, KeyLength, Range);

% Initialize the controller
Sim.Con = Sim.Con.ClearTorques();
Sim.Con.FBType = 2;
Sim.Con.MinSat = [-MaxAnkleT*ones(1,NAnkleT),-MaxHipT*ones(1,NHipT)];
Sim.Con.MaxSat = [MaxAnkleT*ones(1,NAnkleT),MaxHipT*ones(1,NHipT)];

Sim.PMFull = 1;
Sim.Env = Sim.Env.Set('Type','finite','start_slope',start_slope,...
                                      'end_slope',0,...
                                      'parK',0.01);

Sim = Gen.Decode(Sim, Sequence);

% Simulation parameters
% Sim = Sim.SetTime(0,40);
% Sim.RSkip = 150;
Sim = Sim.SetTime(0,0.07,100);

% Set internal parameters (state dimensions, events, etc)
Sim = Sim.Init();

% Some more simulation initialization
Sim.Mod.LegShift = Sim.Mod.Clearance;
Sim.Con = Sim.Con.Reset(Sim.IC(Sim.ConCo));
Sim.Con = Sim.Con.HandleExtFB(Sim.IC(Sim.ModCo),...
            Sim.IC(Sim.ConCo),Sim.Env.SurfSlope(Sim.Mod.xS));

% Simulate
Sim = Sim.Run();

% Calculate eigenvalues
if Sim.Out.Type == 5
    [EigVal,EigVec] = Sim.Poincare();
    % Do some plots
    disp(EigVal);
    disp(abs(EigVal));
    % Display steady state initial conditions
    disp(Sim.ICstore(:,1)');
else
    EigVal = 2*ones(4,1);
    disp(Sim.Out.Text);
end

end