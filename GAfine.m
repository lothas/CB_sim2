function [  ] = GAfine(  )
% Run MOOGA to further fine tune the sequences found in stages
% close all; clear all; clear classes;

% % Prepare input file
% In1 = load('GA_07_23_17_09');
% In2 = load('GA_07_19_00_57');
% GA = In1.GA;
% GA.Seqs = [GA.Seqs(:,:,GA.Progress); In2.GA.Seqs(:,:,In2.GA.Progress)];
% GA.Fit = [GA.Fit(:,:,GA.Progress); In2.GA.Fit(:,:,In2.GA.Progress)];
% GA.Progress = 1;
% save('GA_combined_in','GA');

GA = MOOGA(30,10000);
GA.Fittest = [1900,1900,200];
% GA = MOOGA(50,1500);
% GA.Fittest = [250,250,10];
% GA.FileIn = 'GA_combined_in.mat';
GA.FileIn = 'GA_08_02_22_10.mat';
% GA.FileOut = GA.FileIn;

GA.FileOut = ['GA_',datestr(now,'mm_dd_hh_MM'),'.mat'];
GA.Graphics = 0;

% GA.ReDo = 1;

% Set up the genome
% Controller with push-off, swing pulse + limited ankle pulse
% and feedback
NAnkleT = 1;
NHipT = 5;
MaxAnkleT = 10;
MaxHipT = 50;
MaxTO = 2000; % Max toe off "impulse"
TorqueFBMin = [-2000,-50*ones(1,NAnkleT),-300*ones(1,NHipT)];
TorqueFBMax = -TorqueFBMin;
if NAnkleT>0
    Keys = {'omega0','P_LegE','ExtPulses','Pulses','Pulses',...
        'kOmega_u','kOmega_d','kTorques_u','kTorques_d';...
                   1,       1,      [1,1],   [NAnkleT,1],   [NHipT,2],...
                 1,         1,           1,           1};
    Range = {0.5, 0.55, [-MaxTO, 0.005], [-MaxAnkleT, 0, 0.01], [-MaxHipT, 0, 0.01],...
    -10, -10, TorqueFBMin, TorqueFBMin; % Min
           2, 0.85, [0000, 0.005], [MaxAnkleT, 0.99, 0.99], [MaxHipT, 0.99, 0.99],...
     10,  10, TorqueFBMax, TorqueFBMax}; % Max
else
    Keys = {'omega0','P_LegE','ExtPulses','Pulses',...
    'kOmega_u','kOmega_d','kTorques_u','kTorques_d';...
               1,       1,      [1,1],   [NHipT,2],...
             1,         1,           1,           1};
    Range = {0.5, 0.55, [-MaxTO, 0.005], [-MaxHipT, 0, 0.01],...
    -10, -10, TorqueFBMin, TorqueFBMin; % Min
           2, 0.85, [0000, 0.005], [MaxHipT, 0.99, 0.99],...
     10,  10, TorqueFBMax, TorqueFBMax}; % Max
end
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
GA.Sim.Mod = GA.Sim.Mod.Set('damp',0,'I',0);

% Set up the terrain
start_slope = 0;
GA.Sim.Env = GA.Sim.Env.Set('Type','inc','start_slope',start_slope);

% Initialize the controller
GA.Sim.Con = GA.Sim.Con.ClearTorques();
GA.Sim.Con.FBType = 0;
GA.Sim.Con.MinSat = [-MaxTO,-MaxAnkleT*ones(1,NAnkleT),-MaxHipT*ones(1,NHipT)];
GA.Sim.Con.MaxSat = [0,MaxAnkleT*ones(1,NAnkleT),MaxHipT*ones(1,NHipT)];

% Simulation parameters
GA.Sim.IC = [start_slope, start_slope, 0, 0, 0];
GA.Sim = GA.Sim.SetTime(0,0.15,40);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;
GA.Sim.Con = GA.Sim.Con.HandleEvent(1, GA.Sim.IC(GA.Sim.ConCo));
% GA.Sim.Con = GA.Sim.Con.HandleExtFB(GA.Sim.IC(GA.Sim.ModCo),...
%                                           GA.Sim.IC(GA.Sim.ConCo));
                                
% Fitness functions
GA.FitFcn = {@GA.VelFit;
             @GA.NrgEffFit;
             @GA.EigenFit;
             @GA.UpSlopeFit; % @GA.UphillFitRun;
             @GA.DownSlopeFit; % @GA.DownhillFitRun};
             @GA.ZMPFit};
GA.NFit = length(GA.FitFcn);
GA.Sim.PMFull = 1; % Run poincare map on all 5 coords

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
% GAS = load('TestGA3_07_08_15_53');
% GA.Seqs(end-99:end,:,1) = repmat(GAS.GA.BaseSeq,100,1);
% GA.Seqs(end-99:end,11:18,1) = GAS.GA.Seqs(GAS.GA.GetTopPop(100),:,end);

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
