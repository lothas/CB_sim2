function [  ] = TestGAstages(  )
% Run simulation in stages
% First stage: large population to find solutions faster
% Second stage: small population to evolve faster
% Third/Fourth stage: evolving feedback gains for up/downhill
nStages = 4;
GA = cell(nStages,1);
GAstages = [20, 1000;
            30, 250;
            50, 50;
            50, 250];

time_stamp = datestr(now,'mm_dd_hh_MM');
Graphics = 0;
for s = 1:nStages
    GA{s} = MOOGA(GAstages(s,1),GAstages(s,2));
    GA{s}.FileOut = ['TestGA',int2str(s),'_',time_stamp,'.mat'];
    GA{s}.Graphics = Graphics;
end
GA{1}.FileIn = 'TestGA1_06_25_18_38.mat';
GA{2}.FileIn = 'TestGA2_06_25_18_38.mat';
% GA{3}.FileIn = 'TestGA3_06_23_16_00.mat';

% GA.ReDo = 1;
GenType = 2;

% Set up the genome
switch GenType
    case 1
        % Pulses controller (only 1 pulse to hip)
        Keys = {'omega0','P_LegE','Pulses';...
                       1,       1,   [1,2]};
        Range = {0.5, 0.55, [-20, 0, 0.01]; % Min
                   2, 0.85, [20, 0.99, 0.99]}; % Max
    case 2
        % Controller with push-off, swing pulse + limited ankle pulse
        Keys = {'omega0','P_LegE','ExtPulses','Pulses','Pulses';...
                       1,       1,      [1,1],   [1,1],   [1,2]};
        Range = {0.5, 0.55, [-30, 0.005], [-2, 0, 0.01], [0, 0, 0.01]; % Min
                   2, 0.85, [0, 0.005], [2, 0.99, 0.99], [20, 0.99, 0.99]}; % Max
        Keys2 = {'kOmega_u','kTorques_u';...
                          1,           1};
        Range2 = {-200, [-200,-50,-200]; % Min
                   200, [200,50,200]}; % Max
        Keys3 = {'kOmega_d','kTorques_d';...
                          1,           1};
end
GA{1}.Gen = Genome(Keys, Range);
GA{2}.Gen = Genome(Keys, Range);
KeyLength = GA{2}.Gen.KeyLength;
KeyLength.kTorques_u = 3;
KeyLength.kTorques_d = 3;
GA{3}.Gen = Genome(Keys2, KeyLength, Range2);
GA{4}.Gen = Genome(Keys3, KeyLength, Range2);

for s = 1:nStages
    % Set up the simulations
    GA{s}.Sim = Simulation();
    GA{s}.Sim.Graphics = GA{s}.Graphics;
    GA{s}.Sim.EndCond = 2; % Run until converge (or fall)
    
    % Set up the compass biped model
    GA{s}.Sim.Mod = GA{s}.Sim.Mod.Set('damp',0,'I',0);
end

% Set up the terrain
start_slope = 0;
for s = 1:nStages-2
    GA{s}.Sim.Env = GA{s}.Sim.Env.Set('Type','inc','start_slope',start_slope);
end

% Initialize the controller
for s = 1:nStages
    GA{s}.Sim.Con = GA{s}.Sim.Con.ClearTorques();
    GA{s}.Sim.Con.FBType = 0;
    GA{s}.Sim.IC = [start_slope, start_slope, 0, 0, 0];

    % Simulation parameters
    GA{s}.Sim = GA{s}.Sim.SetTime(0,0.15,30);

    % Some more simulation initialization
    GA{s}.Sim.Mod.LegShift = GA{s}.Sim.Mod.Clearance;
    % Sim.Con = Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));
    GA{s}.Sim.Con = GA{s}.Sim.Con.HandleExtFB(GA{s}.Sim.IC(GA{s}.Sim.ModCo),...
                                              GA{s}.Sim.IC(GA{s}.Sim.ConCo));
end
                                
% Fitness functions
for s = 1:nStages-2
    ThisGA = GA{s};
    ThisGA.NFit = 4;
    ThisGA.FitFcn = {@ThisGA.HeightFit;
                     @ThisGA.VelFit;
                     @ThisGA.NrgEffFit;
                     @ThisGA.EigenFit};
    GA{s} = ThisGA;
end

ThisGA = GA{end-1};
ThisGA.NFit = 1;
ThisGA.FitFcn = {@ThisGA.UphillFit};
ThisGA.Sim.Con.FBType = 1;
GA{end-1} = ThisGA;
ThisGA = GA{end};
ThisGA.NFit = 1;
ThisGA.FitFcn = {@ThisGA.DownhillFit};
ThisGA.Sim.Con.FBType = 1;
GA{end} = ThisGA;

% Run first stage
disp('Running stage 1');
GA{1} = GA{1}.InitGen();
GA{1} = GA{1}.Run();
GA{1}.Plot('Fit');

% Run second stage
disp('Running stage 2');
if isempty(GA{2}.FileIn)
    if exist(GA{1}.FileOut,'file')
        GA{2}.FileIn = GA{1}.FileOut;
    else
        GA{2}.FileIn = GA{1}.FileIn;
    end
    GA{2}.ReDo = 1;
end
GA{2} = GA{2}.InitGen();
GA{2} = GA{2}.Run();
GA{2}.Plot('Fit');

% Select best genomes for stage 3 and 4
Fitness = GA{nStages-2}.Fit(:,:,GA{nStages-2}.Progress);
Weights = [0.0; % height
           0.3; % velocity
           0.2; % energy
           0.5]; % stability
WeiFit = [(1:GA{nStages-2}.Population)',Fitness*Weights];
SortedFit = sortrows(WeiFit,-2);
% Select the best genome from SortedFit
BestID = SortedFit(1,1);

for s = nStages-1:nStages
    GA{s}.BaseGen = GA{nStages-2}.Gen;
    GA{s}.BaseSeq = GA{nStages-2}.Seqs(BestID,:,GA{nStages-2}.Progress);
    GA{s}.BaseSeq(2) = 0.55;
end

% Set the slope in the simulations for the up/downhill stages
leadway = min(max(0.8/Fitness(BestID,2),2),5);
parK = min(max(0.025*Fitness(BestID,2),0.005),0.03);
GA{end-1}.Sim.Env = ...
    GA{end-1}.Sim.Env.Set('Type','inf','start_slope',0,'parK',parK,'start_x',leadway);
GA{end}.Sim.Env = ...
    GA{end}.Sim.Env.Set('Type','inf','start_slope',0,'parK',-parK,'start_x',leadway);

% Run third and fourth stages
disp('Running stage 3');
GA{3} = GA{3}.InitGen();
GA{3} = GA{3}.Run();
GA{3}.Plot('Fit');

disp('Running stage 4');
GA{4} = GA{4}.InitGen();
GA{4} = GA{4}.Run();
GA{4}.Plot('Fit');



