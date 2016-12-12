function [  ] = GA_Vel_Ours( gen, pop, file_in, file_out )
% Run MOOGA using only Vel, Nrg and Speed-range fitness but limiting the
% simulation with a bounded foot size (ZMP threshold)

if nargin<4
    GA = MOOGA(30,750);
%     GA = MOOGA(10,500);
    GA = GA.SetFittest(15,15,0.5);
    GA.JOAT = 0; GA.Quant = 0.8;
    % GA.Fittest = [20,20,1];
%     GA.FileIn = 'VGAM_11_22_00_00.mat';
%     GA.FileOut = GA.FileIn;

    GA.FileOut = ['VGAO_',datestr(now,'mm_dd_hh_MM'),'.mat'];
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
genome_file = 'SpitzGenome.mat';
if exist(genome_file, 'file') ~= 2
    nAnkle = 2; % Number of ankle torques
    nHip = 3;   % Number of hip torques
    maxAnkle = 8;   % Max ankle torque
    maxHip = 20;    % Max hip torque
    N = nAnkle+nHip; %#ok<NASGU>
    
    TorqueFBMin = [-3*ones(1,nAnkle),-6*ones(1,nHip)];
    TorqueFBMax = -TorqueFBMin;
    Keys = {'Freq',    'Pulses',   'Pulses',  'kFreq','kTorques';...
                 1,   [nAnkle,1],  [nHip,2],        1,         1};
    Range = {  1.2, [-maxAnkle, 0, 0.01], [-maxHip, 0, 0.01],...
                                                  -1, TorqueFBMin; % Min
               1.6, [maxAnkle, 0.99, 0.99], [maxHip, 0.99, 0.99],...
                                                   1, TorqueFBMax}; % Max

    MutDelta0 = 0.04;
    MutDelta1 = 0.02;
    
    save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
        'N', 'TorqueFBMin', 'TorqueFBMax', ...
        'MutDelta0', 'MutDelta1', 'Keys', 'Range');
else
    load(genome_file);
end

GA.Gen = Genome(Keys, Range);
KeyLength = GA.Gen.KeyLength;
KeyLength.kTorques = length(TorqueFBMin);
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
GA.Sim.Con.MinSat = [-maxAnkle*ones(1,nAnkle),-maxHip*ones(1,nHip)];
GA.Sim.Con.MaxSat = [maxAnkle*ones(1,nAnkle),maxHip*ones(1,nHip)];

% Simulation parameters
GA.Sim.IC = [start_slope, start_slope, 0, 0, 0];
GA.Sim = GA.Sim.SetTime(0,0.03,20);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;
% GA.Sim.Con = GA.Sim.Con.HandleEvent(1, GA.Sim.IC(GA.Sim.ConCo));
% GA.Sim.Con = GA.Sim.Con.HandleExtFB(GA.Sim.IC(GA.Sim.ModCo),...
%                                           GA.Sim.IC(GA.Sim.ConCo));
                                
% Fitness functions
GA.FitFcn = {1, @MOOGA.VelFit;
             2, @MOOGA.NrgEffFit;
             3:10, @MOOGA.VelRangeFit;
             11, @MOOGA.EigenFit};
GA.FitIDs = [1,2,3]; % Velocity and average COT
GA.FitMinMax = [1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1];
% Pareto always looks to maximize the fitness value but sometimes we want
% to check values that don't go into the pareto selection and we want to
% minimize (or get the maximum negative value, e.g. a negative slope)

GA.NFit = size(GA.FitIDs,2);
GA.Sim.PMFull = 1; % Run poincare map on all coords

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
GA.PrintRuntime();
end
