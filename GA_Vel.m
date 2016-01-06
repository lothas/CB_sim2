function [  ] = GA_Vel( gen, pop, file_in, file_out )
% Run MOOGA using only Vel, Nrg and Speed-range fitness but limiting the
% simulation with a bounded foot size (ZMP threshold)

if nargin<4
    GA = MOOGA(30,1200);
    % GA = MOOGA(10,100);
    GA = GA.SetFittest(15,15,0.5);
    GA.JOAT = 2; GA.Quant = 0.5;
    % GA.Fittest = [20,20,1];
%     GA.FileIn = 'VGA_12_14_08_57.mat';
%     GA.FileOut = GA.FileIn;

    GA.FileOut = ['VGA_',datestr(now,'mm_dd_hh_MM'),'.mat'];
else
    GA = MOOGA(gen,pop);
    GA = GA.SetFittest(20,20,0.5);
    GA.JOAT = 2; GA.Quant = 0.6;
    % GA.Fittest = [20,20,1];
    GA.FileIn = file_in;
    GA.FileOut = file_out;
end

GA.Graphics = 0;
GA.ReDo = 1;

% Set up the genome
% Controller with swing pulse + limited ankle pulse and feedback
NAnkleT = 1;
NHipT = 2;
MaxAnkleT = 25;
MaxHipT = 25;
TorqueFBMin = [-5*ones(1,NAnkleT),-10*ones(1,NHipT)];
TorqueFBMax = -TorqueFBMin;
Keys = {'omega0',    'Pulses',   'Pulses',...
      'sOmega_f','sTorques_f', 'sOmega_s','sTorques_s';...
               1, [NAnkleT,1],  [NHipT,2],...
               1,           1,          1,          1};
Range = {0.7, [-MaxAnkleT, 0, 0.01], [-MaxHipT, 0, 0.01],...
    -0.3, TorqueFBMin, -6, TorqueFBMin; % Min
         1.7, [MaxAnkleT, 0.99, 0.99], [MaxHipT, 0.99, 0.99],...
    0.3, TorqueFBMax,  6, TorqueFBMax}; % Max
MutDelta0 = 0.04;
MutDelta1 = 0.02;

GA.Gen = Genome(Keys, Range);
KeyLength = GA.Gen.KeyLength;
KeyLength.sTorques_f = length(TorqueFBMax);
KeyLength.sTorques_s = length(TorqueFBMin);
GA.Gen = Genome(Keys, KeyLength, Range);

% Set up the simulations
GA.Sim = Simulation();
GA.Sim.Graphics = GA.Graphics;
GA.Sim.EndCond = 2; % Run until converge (or fall)

% Set up the compass biped model
% GA.Sim.Mod = GA.Sim.Mod.Set('damp',0.3);
GA.Sim.Mod = GA.Sim.Mod.Set('I',0,'damp',0,'A2T',0.16,'A2H',0.12);

% Set up the terrain
start_slope = 0;
GA.Sim.Env = GA.Sim.Env.Set('Type','inc','start_slope',start_slope);

% Initialize the controller
GA.Sim.Con = GA.Sim.Con.ClearTorques();
GA.Sim.Con.FBType = 0; % no slope feedback
GA.Sim.Con.MinSat = [-MaxAnkleT*ones(1,NAnkleT),-MaxHipT*ones(1,NHipT)];
GA.Sim.Con.MaxSat = [MaxAnkleT*ones(1,NAnkleT),MaxHipT*ones(1,NHipT)];

% Simulation parameters
GA.Sim.IC = [start_slope, start_slope, 0, 0, 0];
GA.Sim = GA.Sim.SetTime(0,0.15,40);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;
% GA.Sim.Con = GA.Sim.Con.HandleEvent(1, GA.Sim.IC(GA.Sim.ConCo));
% GA.Sim.Con = GA.Sim.Con.HandleExtFB(GA.Sim.IC(GA.Sim.ModCo),...
%                                           GA.Sim.IC(GA.Sim.ConCo));
                                
% Fitness functions
GA.FitFcn = {1, @MOOGA.VelFit;
             2, @MOOGA.NrgEffFit;
             3:7, @MOOGA.VelRangeFit;
             8, @MOOGA.EigenFit};
GA.FitIDs = [1,2,3];
GA.NFit = size(GA.FitIDs,2);
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
% TrySim = GA.Gen.Decode(TrySim,GA.Seqs(1,:,1));
% TrySim = TrySim.SetTime(0,0.01,10);
% TrySim = TrySim.Init();
% TrySim.Run()

GA = GA.Run();
GA.Plot('Fit');

end
