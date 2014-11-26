function [  ] = GA_SlopeFB(  )
% Run MOOGA to fine tune the hand tuned controller used in the ROBIO/TRO
% paper.
% close all; clear all; clear classes;

GA = MOOGA(15,1000);
% GA = MOOGA(10,100);
GA = GA.SetFittest(20,20,0.5);
GA.JOAT = 2; GA.Quant = 0.6;
% GA.Fittest = [20,20,1];
% GA.FileIn = 'GA_09_05_11_08.mat';
% GA.FileOut = GA.FileIn;

GA.FileOut = ['GA_',datestr(now,'mm_dd_hh_MM'),'.mat'];
GA.Graphics = 1;

% GA.ReDo = 1;

% Set up the base genome
NAnkleT = 3;
NHipT = 3;
MaxAnkleT = 150;
MaxHipT = 150;
Keys = {'omega0','P_LegE',    'Pulses',    'Pulses';...
               1,       1, [NAnkleT,1],   [NHipT,2]};
Range = {0.7, 0.55, [-MaxAnkleT, 0, 0.01], [-MaxHipT, 0, 0.01]; % Min
         1.5, 0.85, [MaxAnkleT, 0.99, 0.99], [MaxHipT, 0.99, 0.99]}; % Max
GA.BaseGen = Genome(Keys, Range);
GA.BaseSeq = [1.47469, 0.727897,...
            -24.2109, 0.539044, 0.201354,... % Ankle pulse 1
            8.91417, 0.532459, 0.259763,... % Ankle pulse 2
            4.85994, 0.591071, 0.15863,... % Ankle pulse 3
            56.5475, 0.438046, 0.193643,... % Hip pulse 1
           -50.3699, 0.462751, 0.17404,... % Hip pulse 2
           -13.8021, 0.560096, 0.297654]; % Hip pulse 3

% Set up the genome
% Controller with feedback
TorqueFBMin = [-1000*ones(1,NAnkleT),-1000*ones(1,NHipT)];
TorqueFBMax = -TorqueFBMin;
Keys = {'kOmega_u','kOmega_d','kTorques_u','kTorques_d';
                1,         1,           1,           1};
Range = {0, 0, TorqueFBMin, TorqueFBMin; % Min
         3, 3, TorqueFBMax, TorqueFBMax}; % Max

GA.Gen = Genome(Keys, Range);
KeyLength = GA.Gen.KeyLength;
KeyLength.kTorques_u = length(TorqueFBMin);
KeyLength.kTorques_d = length(TorqueFBMin);
GA.Gen = Genome(Keys, KeyLength, Range);
MutDelta0 = 0.03;
MutDelta1 = 0.01;

% Set up the simulations
GA.Sim = Simulation();
GA.Sim.Graphics = GA.Graphics;
GA.Sim.EndCond = 2; % Run until converge (or fall)

% Set up the compass biped model
GA.Sim.Mod = GA.Sim.Mod.Set('damp',0,'I',0);

% Set up the terrain
start_slope = 0;
GA.Sim.Env = GA.Sim.Env.Set('Type','inc','start_slope',start_slope);

% Initialize the controller
GA.Sim.Con = GA.Sim.Con.ClearTorques();
GA.Sim.Con.FBType = 0;
GA.Sim.Con.MinSat = [-MaxAnkleT*ones(1,NAnkleT),-MaxHipT*ones(1,NHipT)];
GA.Sim.Con.MaxSat = [MaxAnkleT*ones(1,NAnkleT),MaxHipT*ones(1,NHipT)];

% Simulation parameters
GA.Sim.IC = [start_slope, start_slope, 0, 0, 0];
GA.Sim = GA.Sim.SetTime(0,0.15,10);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;
% GA.Sim.Con = GA.Sim.Con.HandleEvent(1, GA.Sim.IC(GA.Sim.ConCo));
% GA.Sim.Con = GA.Sim.Con.HandleExtFB(GA.Sim.IC(GA.Sim.ModCo),...
%                                           GA.Sim.IC(GA.Sim.ConCo));
                                
% Fitness functions
GA.FitFcn = {1:3, @MOOGA.ZMPUpFit;
             4:6, @MOOGA.ZMPDownFit};
GA.FitIDs = [1,4];
GA.NFit = size(GA.FitFcn,1);
GA.Sim.PMFull = 1; % Run poincare map on all 5 coords

GA = GA.InitGen();

% Update MOOGA parameters after each generation
    function GA = GenFcn(GA)
        % Increase allowed genome range
        % Interpolate range
%         i = min(GA.Progress,10)/10;
%         Range = i*Range2 + (1-i)*Range1;
%         GA.Gen.Range = Range;
        
        % Reduce mutation range
%         j = min(GA.Progress,30)/30;
        j = GA.Progress/GA.Generations;
        GA.Gen.MutDelta = (1-j)*MutDelta0 + MutDelta1*j;
    end
GA.GenerationFcn = @GenFcn;

% TrySim = deepcopy(GA.Sim);
% TrySim = GA.BaseGen.Decode(TrySim,GA.BaseSeq);
% TrySim = GA.Gen.Decode(TrySim,GA.Seqs(1,:,1));
% TrySim = TrySim.SetTime(0,0.01,10);
% TrySim = TrySim.Init();
% TrySim.Run()

GA = GA.Run();
GA.Plot('Fit');

end
