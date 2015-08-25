function [  ] = GA_Slope( gen, pop, slope, file_in, file_out )
% Run MOOGA using Vel, Nrg, Robustness and Uniqueness fitness
% The Simulation is limited with a bounded foot size (ZMP threshold)
% Inputs:
% gen - number of generations to run
% pop - number of genomes in population
% slope - desired slope for optimization (in degrees)
% file_in - input file with initial population / to continue MOGA
% file_out - output file for results

if nargin<5
    GA = MOOGA(25,1000);
    % GA = MOOGA(10,100);
    GA = GA.SetFittest(15,15,0.5);
    slope = 1;
    % GA.Fittest = [20,20,1];
    GA.FileIn = 'GA_0_08_22_17_51.mat';
%     GA.FileOut = GA.FileIn;

    GA.FileOut = ['GA_',num2str(slope),'_',datestr(now,'mm_dd_hh_MM'),'.mat'];
else
    GA = MOOGA(gen,pop);
    GA = GA.SetFittest(15,15,0.5);
    % GA.Fittest = [20,20,1];
    GA.FileIn = file_in;
    GA.FileOut = file_out;
end

GA.Graphics = 0;
GA.ReDo = 1; % 1 - start the algorithm from scratch (random pop)

% Set up the genome
% Controller with swing pulse + limited ankle pulse and feedback
NAnkleT = 1;
NHipT = 2;
MaxAnkleT = 20;
MaxHipT = 80;
Keys = {'IC', 'omega0',    'Pulses',    'Pulses';
           4,        1, [NAnkleT,1],   [NHipT,2]};
Range = {[0,-2,-2,0], 0.90, [-MaxAnkleT, 0, 0.01], [-MaxHipT, 0, 0.01];... % Min
     [0.79,2,2,0.99], 1.8, [MaxAnkleT, 0.99, 0.99], [MaxHipT, 0.99, 0.99]}; % Max
MutDelta0 = 0.04;
MutDelta1 = 0.02;

GA.Gen = Genome(Keys, Range);

% Set up the simulations
GA.Sim = Simulation();
GA.Sim.Graphics = GA.Graphics;
GA.Sim.EndCond = 2; % Run until converge (or fall)

% Set up the compass biped model
% GA.Sim.Mod = GA.Sim.Mod.Set('damp',0.3);
GA.Sim.Mod = GA.Sim.Mod.Set('I',0,'damp',0,'A2T',0.16,'A2H',0.12);

% Set up the terrain
GA.Sim.Env = Terrain(0,slope);

% Initialize the controller
GA.Sim.Con = GA.Sim.Con.ClearTorques();
GA.Sim.Con.FBType = 0;  % no feedback
GA.Sim.Con.MinSat = [-MaxAnkleT*ones(1,NAnkleT),-MaxHipT*ones(1,NHipT)];
GA.Sim.Con.MaxSat = [MaxAnkleT*ones(1,NAnkleT),MaxHipT*ones(1,NHipT)];

% Simulation parameters
GA.Sim.IC = [slope*pi/180, slope*pi/180, 0, 0, 0];
tstep = 0.05;  % simulation time step
tend = 40;  % max simulation runtime
GA.Sim = GA.Sim.SetTime(0,tstep,tend);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;
                                
% Fitness functions
% GA.FitFcn = {1, @MOOGA.VelFit;
%              2, @MOOGA.NrgEffFit;
%              3, @MOOGA.RobustFit;
%              4, @MOOGA.EigenFit};
GA.FitFcn = {1, @MOOGA.VelFit;
             2, @MOOGA.NrgEffFit;
             3, @MOOGA.EigenFit};
GA.FitIDs = [1,2,3];
GA.NFit = size(GA.FitFcn,1);
GA.Sim.PMFull = 1; % Run poincare map on all 5 coords
GA.SelectionMode = 3;  % Run "crowd control" to prevent solutions from
                       % clustering

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
