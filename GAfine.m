function [  ] = GAfine(  )
% Run MOOGA to further fine tune the sequences found in stages
% close all; clear all; clear classes;

GA = MOOGA(25,2500);
GA.FileIn = 'GA_07_04_04_32.mat';
% GA.FileIn = 'GA_07_09_16_33.mat';
% GA.FileOut = GA.FileIn;
GA.FileOut = ['GA_',datestr(now,'mm_dd_hh_MM'),'.mat'];
GA.Graphics = 0;

GA.ReDo = 1;

% Set up the genome
% Controller with push-off, swing pulse + limited ankle pulse
% and feedback
Keys = {'omega0','P_LegE','ExtPulses','Pulses','Pulses',...
    'kOmega_u','kOmega_d','kTorques_u','kTorques_d';...
               1,       1,      [1,1],   [1,1],   [1,2],...
             1,         1,           1,           1};
Range = {0.5, 0.55, [-200, 0.005], [-40, 0, 0.01], [-200, 0, 0.01],...
    -200, -200, [-400,-200,-400], [-400,-200,-400]; % Min
           2, 0.85, [200, 0.005], [40, 0.99, 0.99], [200, 0.99, 0.99],...
     200,  200, [400,200,400], [400,200,400]}; % Max
Range1 = {0.5, 0.55, [-30, 0.005], [-2, 0, 0.01], [-20, 0, 0.01],...
    -200, -200, [-200,-50,-200], [-200,-50,-200]; % Min
           2, 0.85, [30, 0.005], [2, 0.99, 0.99], [20, 0.99, 0.99],...
     200,  200, [200,50,200], [200,50,200]}; % Max
Range2 = {0.5, 0.55, [-200, 0.005], [-40, 0, 0.01], [-200, 0, 0.01],...
    -200, -200, [-400,-200,-400], [-400,-200,-400]; % Min
           2, 0.85, [200, 0.005], [40, 0.99, 0.99], [200, 0.99, 0.99],...
     200,  200, [400,200,400], [400,200,400]}; % Max
Range3 = {0.5, 0.55, [-400, 0.005], [-100, 0, 0.01], [-400, 0, 0.01],...
    -200, -200, [-800,-400,-800], [-800,-400,-800]; % Min
           2, 0.85, [400, 0.005], [100, 0.99, 0.99], [400, 0.99, 0.99],...
     200,  200, [800,400,800], [800,400,800]}; % Max
MutDelta0 = 0.1;
MutDelta1 = 0.01;

GA.Gen = Genome(Keys, Range3);
KeyLength = GA.Gen.KeyLength;
KeyLength.kTorques_u = 3;
KeyLength.kTorques_d = 3;
Range2 = GA.Gen.Range;
GA.Gen = Genome(Keys, KeyLength, Range3);
Range1 = GA.Gen.Range;

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
GA.Sim.IC = [start_slope, start_slope, 0, 0, 0];

% Simulation parameters
GA.Sim = GA.Sim.SetTime(0,0.15,40);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;
% Sim.Con = Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));
GA.Sim.Con = GA.Sim.Con.HandleExtFB(GA.Sim.IC(GA.Sim.ModCo),...
                                          GA.Sim.IC(GA.Sim.ConCo));
GA.Sim = GA.Sim.Init();
                                
% Fitness functions
GA.NFit = 6;
GA.FitFcn = {@GA.VelFit;
             @GA.NrgEffFit;
             @GA.EigenFit;
             @GA.UphillFitRun;
             @GA.DownhillFitRun;
             @GA.ZMPFit};

GA = GA.InitGen();
% Add the best guesses from previous runs
% GA.Seqs(1:10,1:10,1) = ...
%     [1.32491,0.66107,-0.01570,0.005,-1.98446,0.01223,0.76611,7.56098,0.06446,0.10987;
%      1.32491,0.72147,-0.01570,0.005,-1.98446,0.01223,0.76611,7.56098,0.06446,0.10987;
%      1.32491,0.72147,-0.01570,0.005,-1.98446,0.01223,0.76611,7.75285,0.06446,0.10987;
%      1.33672,0.74226,-0.28342,0.005,-1.78457,0.12747,0.76611,7.19713,0.09427,0.10806;
%      1.34743,0.65412,-0.01570,0.005,-1.75601,0.07494,0.72612,5.71611,0.00186,0.11990;
%      1.33672,0.74226,-0.28342,0.005,-1.78457,0.08475,0.76611,7.80237,0.07867,0.10806;
%      1.34743,0.68067,-0.01570,0.005,-1.75601,0.03990,0.72612,4.72148,0.00186,0.15320;
%      1.33078,0.69347,-0.01570,0.005,-1.98446,0.01223,0.76611,4.72148,0.00186,0.15552;
%      1.32491,0.63032,-0.01570,0.005,-1.87480,0.04302,0.78613,10.25600,0.06447,0.06590;
%      1.31369,0.72147,-0.01570,0.005,-1.98446,0.01727,0.81724,6.20295,0.11312,0.10987];
%  
% DownGains = ...
%     [161.30011,-48.56199,73.562426;
%      197.29988,-49.59438,73.962086;
%      197.29988,-49.59438,73.562426;
%      193.57764,-49.59438,73.562426;
%      193.57764,-49.59438,73.562426;
%      193.57764,-49.59438,73.562426;
%      192.85560,-49.59438,73.583720;
%      193.57764,-49.59438,73.562426;
%      193.57764,-49.59438,73.562426;
%      178.82647,-49.28193,74.315666];
%  
% GA.Seqs(1:10,11:18,1) = [repmat(1.09517,10,1),...
%                          repmat(2.43205,10,1),...
%                          repmat([-163.33421,-47.58048,154.71160],10,1),...
%                          DownGains];

% Also add the best genes from TestGAStages
GAS = load('TestGA3_07_08_15_53');
GA.Seqs(end-99:end,:,1) = repmat(GAS.GA.BaseSeq,100,1);
GA.Seqs(end-99:end,11:18,1) = GAS.GA.Seqs(GAS.GA.GetTopPop(100),:,end);

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

GA = GA.Run();
GA.Plot('Fit');

end
