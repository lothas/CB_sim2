function [  ] = GA_Vel_Matsuoka( gen, pop, file_in, file_out )
% Run MOOGA using only Vel, Nrg and Speed-range fitness but limiting the
% simulation with a bounded foot size (ZMP threshold)

if nargin<4
    GA = MOOGA(10,6000);
%     GA = MOOGA(10,100);
    GA = GA.SetFittest(15,15,0.5);
    GA.JOAT = 2; GA.Quant = 0.5;
    % GA.Fittest = [20,20,1];
    GA.FileIn = 'VGAM_02_01_21_36.mat';
%     GA.FileOut = GA.FileIn;

    GA.FileOut = ['VGAM_',datestr(now,'mm_dd_hh_MM'),'.mat'];
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
nAnkle = 1; % Number of ankle torques
nHip = 1;   % Number of hip torques
maxAnkle = 8;   % Max ankle torque
maxHip = 20;    % Max hip torque
Mamp = [maxAnkle*ones(1,2*nAnkle), maxHip*ones(1,2*nHip)];
mamp = 0*Mamp;
N = nAnkle+nHip;
Mwex = -10*ones(1,4*N);
mwex = 0*Mwex;

% TorqueFBMin = [-5*ones(1,NAnkleT),-10*ones(1,NHipT)];
% TorqueFBMax = -TorqueFBMin;

% Keys = {'tau','tav','beta','amp','win','wex','ks_tau',  'ks_out','IC_matsuoka';
%            1 ,   1 ,    1 , 2*N ,   1 , 2*N ,      1 ,      2*N ,           0 };
% Range = {0.1 , 0.1 ,   10 , mamp, -20 , Mwex,   -1e2 , -0.1*Mamp; % Min
%            2 ,   2 ,   30 , Mamp,  -1 , mwex,    1e2 ,  0.1*Mamp}; % Max
Keys = {'tau','tav','beta','amp','win','wex','ks_tau',  'ks_out','IC_matsuoka';
           1 ,   1 ,    1 , 2*N ,   1 , 4*N ,      1 ,      2*N ,           0 };
Range = {0.1 , 0.1 ,    1 , mamp, -10 , Mwex,   -1e2 , -0.1*Mamp; % Min
           2 ,   2 ,   30 , Mamp, -.5 , mwex,    1e2 ,  0.1*Mamp}; % Max

MutDelta0 = 0.04;
MutDelta1 = 0.02;

GA.Gen = Genome(Keys, Range);
% KeyLength = GA.Gen.KeyLength;
% KeyLength.sTorques_f = length(TorqueFBMax);
% KeyLength.sTorques_s = length(TorqueFBMin);
% GA.Gen = Genome(Keys, KeyLength, Range);

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
GA.Sim.Con = Matsuoka;
GA.Sim.Con.startup_t = 1.5; % Give some time for the neurons to converge
% before applying a torque
GA.Sim.Con.FBType = 0; % no slope feedback
GA.Sim.Con.nPulses = N;
GA.Sim.Con.stDim = 4*N;
GA.Sim.Con = GA.Sim.Con.SetOutMatrix([nAnkle,nHip]);
GA.Sim.Con.MinSat = [-maxAnkle,-maxHip];
GA.Sim.Con.MaxSat = [ maxAnkle, maxHip];

% Simulation parameters
GA.Sim.IC = [start_slope, start_slope, 0, 0, zeros(1, GA.Sim.Con.stDim)];
GA.Sim = GA.Sim.SetTime(0,0.09,40);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;
% GA.Sim.Con = GA.Sim.Con.HandleEvent(1, GA.Sim.IC(GA.Sim.ConCo));
% GA.Sim.Con = GA.Sim.Con.HandleExtFB(GA.Sim.IC(GA.Sim.ModCo),...
%                                           GA.Sim.IC(GA.Sim.ConCo));
                                
% Fitness functions
GA.FitFcn = {1, @MOOGA.VelFit;
             2, @MOOGA.NrgEffFit;
             3, @MOOGA.EigenFit};
GA.FitIDs = [1,2]; % Velocity range and average COT
GA.FitMinMax = [1, 1, 1];

% GA.FitFcn = {1, @MOOGA.VelFit;
%              2, @MOOGA.NrgEffFit;
%              3:7, @MOOGA.VelRangeFit;
%              8, @MOOGA.EigenFit};

% GA.FitFcn = {1, @MOOGA.VelFit;
%              2, @MOOGA.EigenFit;
%              3:10, @MOOGA.VelRangeFit};
% GA.FitIDs = [3,4]; % Velocity range and average COT
% GA.FitMinMax = [1, 1, 1, 1, -1, 1, -1, 1, -1, 1];
% Pareto always looks to maximize the fitness value but sometimes we want
% to check values that don't go into the pareto selection and we want to
% minimize (or get the maximum negative value, e.g. a negative slope)

GA.NFit = size(GA.FitIDs,2);
GA.Sim.PMFull = 1; % Run poincare map on all coords

GA = GA.InitGen();
% GA.Seqs(1,:,1) = [0.3,  0.6,   15, [1.00,4.00,10.00,1.00],     -1,...
%             [-1, 0, 0,-3,-1, 0, 0,-3], -0.0001, [-0.02 0.02 0.4 0.1]];

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
