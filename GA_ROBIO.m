function [  ] = GA_ROBIO(  )
% Run MOOGA to fine tune the hand tuned controller used in the ROBIO/TRO
% paper.
% close all; clear all; clear classes;

GA = MOOGA(35,1600);
% GA = MOOGA(10,100);
GA.Fittest = [300,300,10];
GA.JOAT = 2; GA.Quant = 0.6;
% GA.Fittest = [20,20,1];
% GA.FileIn = 'GA_08_04_15_34.mat';
% GA.FileOut = GA.FileIn;

GA.FileOut = ['GA_',datestr(now,'mm_dd_hh_MM'),'.mat'];
GA.Graphics = 1;

% GA.ReDo = 1;

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
MutDelta0 = 0.1;
MutDelta1 = 0.01;

GA.Gen = Genome(Keys, Range);
KeyLength = GA.Gen.KeyLength;
KeyLength.kTorques_u = length(TorqueFBMin);
KeyLength.kTorques_d = length(TorqueFBMin);
GA.Gen = Genome(Keys, KeyLength, Range);

% Set up the simulations
GA.Sim = Simulation();
GA.Sim.Graphics = GA.Graphics;
GA.Sim.EndCond = 2; % Run until converge (or fall)

% Set up the compass biped model
GA.Sim.Mod = GA.Sim.Mod.Set('damp',0.3);

% Set up the terrain
start_slope = 0;
GA.Sim.Env = GA.Sim.Env.Set('Type','inc','start_slope',start_slope);

% Initialize the controller
GA.Sim.Con = GA.Sim.Con.ClearTorques();
GA.Sim.Con.FBType = 0;
GA.Sim.Con.MinSat = [-MaxAnkleT*ones(1,NAnkleT),-MaxHipT*ones(1,NHipT)];
GA.Sim.Con.MaxSat = [MaxAnkleT*ones(1,NAnkleT),MaxHipT*ones(1,NHipT)];

% Simulation parameters
GA.Sim.IC = [start_slope, start_slope, 0, 0, 0.8];
GA.Sim = GA.Sim.SetTime(0,0.15,40);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;
% GA.Sim.Con = GA.Sim.Con.HandleEvent(1, GA.Sim.IC(GA.Sim.ConCo));
% GA.Sim.Con = GA.Sim.Con.HandleExtFB(GA.Sim.IC(GA.Sim.ModCo),...
%                                           GA.Sim.IC(GA.Sim.ConCo));
                                
% Fitness functions
GA.FitFcn = {@GA.VelFit;
             @GA.NrgEffFit;
             @GA.EigenFit;
             @GA.UpSlopeFit; % @GA.UphillFitRun;
             @GA.DownSlopeFit; % @GA.DownhillFitRun};
             @GA.ZMPFit};
% GA.FitFcn = {@MOOGA.UpSlopeFit; % @GA.UphillFitRun;
%              @MOOGA.DownSlopeFit; % @GA.DownhillFitRun};
%              @MOOGA.ZMPFit};
GA.NFit = length(GA.FitFcn);
GA.Sim.PMFull = 1; % Run poincare map on all 5 coords

GA = GA.InitGen();

% Add the hand-tuned sequence
GA.Seqs(1,:,1) = [1.106, 0.65, 13.626, 0, 0.1, -4.628, 0, 0.4, 0, 0, 0.4,...
                               13.626, 0, 0.1, 0, 0, 0.1, 0, 0, 0.1,...
       0.46, 0.9, 95, -443, 0, 95, 0, 0, 80, -295, 0, 80, 0, 0];

% Update MOOGA parameters after each generation
    function GA = GenFcn(GA)
        % Increase allowed genome range
        % Interpolate range
%         i = min(GA.Progress,10)/10;
%         Range = i*Range2 + (1-i)*Range1;
%         GA.Gen.Range = Range;
        
        % Reduce mutation range
        j = min(GA.Progress,30)/30;
        GA.Gen.MutDelta = (1-j)*MutDelta0 + MutDelta1*j;
    end
GA.GenerationFcn = @GenFcn;

% TrySim = deepcopy(GA.Sim);
% TrySim = GA.Gen.Decode(TrySim,GA.Seqs(1,:,1));
% TrySim = TrySim.SetTime(0,0.01,10);
% TrySim = TrySim.Init();
% TrySim.Run()

GA = GA.Run();
GA.Plot('Fit');

end
