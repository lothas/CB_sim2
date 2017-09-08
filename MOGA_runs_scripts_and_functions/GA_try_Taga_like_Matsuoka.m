function [  ] = GA_try_Taga_like_Matsuoka(whichCase,fileIn)
% Run MOOGA using only Vel, Nrg and Speed-range fitness but limiting the
% simulation with a bounded foot size (ZMP threshold)



% GA = MOOGA(20,500);
GA = MOOGA(10,1000);
GA = GA.SetFittest(15,15,0.5);
GA.JOAT = 2; GA.Quant = 0.7;

GA.FileIn = fileIn;

switch whichCase
    case 'GA only'
        % Use NN?
        use_NN = 0;
        % Rescale?
        GA.rescaleFcn = [];
        % out file name:
        GA.FileOut = ['VGAM_Taga_like_',datestr(now,'mm_dd_hh_MM'),...
            '_GA_only','.mat'];
        
    case 'GA + NN'
        use_NN = 1;
        GA.rescaleFcn = [];
        GA.FileOut = ['VGAM_Taga_like_',datestr(now,'mm_dd_hh_MM'),...
            '_NN_only','.mat'];
    case 'GA + rescale'
        use_NN = 0;
        GA.rescaleFcn = @rescaleFcn;
        GA.FileOut = ['VGAM_Taga_like_',datestr(now,'mm_dd_hh_MM'),...
            '_rescale_only','.mat'];
    case 'GA + NN + rescale'
        use_NN = 1;
        GA.rescaleFcn = @rescaleFcn;
        GA.FileOut = ['VGAM_Taga_like_',datestr(now,'mm_dd_hh_MM'),...
            '_NN_and_rescale','.mat'];
end
        

GA.Graphics = 0;
GA.ReDo = 1;

% Set up the genome
load('MatsuokaGenome_4Neuron_tagaLike.mat','Keys','Range','N',...
    'nAnkle','nHip','maxAnkle', 'maxHip','MutDelta0','MutDelta1');

GA.Gen = Genome(Keys, Range);

MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.nNeurons = 2*N;

% Use NN?
if use_NN
    GANN_file = 'MatsGANN.mat';
    
    maxN = 250000;
    NNSamples = 500;
    
    inFilenames =...
        {'MatsRandomRes_4Neurons_TagaLike_Narrow_b_Narrow_W_All_ocs_1.mat',...
    'MatsRandomRes_4Neurons_TagaLike_Narrow_b_Narrow_All_osc_1_RESCALED.mat'};

    MML.sample_genes = {'\tau_r','4neuron_taga_like'}; % the name of the 'set' options of the Taga like weigths
    MML.target_genes = {'beta'};
    
%     MML.sample_genes = {'4neuron_taga_like'}; 
%     MML.target_genes = {'\tau_r','beta'};

    [samples, targets, normParams] = MML.prepareNNData(inFilenames, maxN);
    MML.normParams = normParams;
    
%     if exist(GANN_file, 'file') ~= 2
        architecture = [20,20];
        [net, ~, ~, ~, ~, ~] = ...
                MML.trainNN(samples, targets, architecture, NNSamples);
        save(GANN_file,'net');
%     else
%         GANN_net = load(GANN_file);
%         net = GANN_net.net;
%     end

    GA.NN = net;
    GA.NNFcn = @NNFcn;
end

function seq = NNFcn(Gen, net, seq, X, T)

    % % don't need to run Sim again. run the NN on all CPG's??
    
    [~, periods, ~, ~, ~] = MML.processResults(X, T);
    % don't do anything if CPG IS stable
    if ~any(isnan(periods)) 
        return
    end

    % Use NN to select best value for tau gene
    desPeriod = MML.perLim(1) + ...
                 rand()*(MML.perLim(2)-MML.perLim(1));

    seq = MML.getNNPar(net, seq, desPeriod);

    ids = seq < Gen.Range(1,:) | seq > Gen.Range(2,:);
    seq(ids) = min(max(seq(ids),Gen.Range(1,ids)),Gen.Range(2,ids));
end

function seq = rescaleFcn(Gen, seq, X, T)
    [~, periods, ~, ~, ~] = MML.processResults(X, T);

    % don't do anything if CPG is not stable
    if any(isnan(periods)) 
        return
    end
    inputPeriod = max(periods);

    % don't do anything if CPG is osc in the right period range
    if (inputPeriod > MML.perLim(1)) && (inputPeriod < MML.perLim(2))
        return
    end

    % Select new random period within desired range
    des_period = MML.perLim(1) + rand()*(MML.perLim(2)-MML.perLim(1));

    % Scale Tr, Ta to obtain desired period
    ratio = des_period/inputPeriod;
    seq(1) = seq(1)*ratio;
    if seq(1) < Gen.Range(1,1) || seq(1) > Gen.Range(2,1)
%             warning('Genetic sequence out of bounds, using bounded tau gene')
        % Bound tau gene
        seq(1) = min(max(seq(1), Gen.Range(1,1)), Gen.Range(2,1));
    end
end

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
GA.Sim.Con.startup_t = 1.0; % Give some time for the neurons to converge
% before applying a torque
GA.Sim.Con.FBType = 0; % no slope feedback
GA.Sim.Con.nPulses = N;
GA.Sim.Con.stDim = 4*N;
GA.Sim.Con = GA.Sim.Con.SetOutMatrix([nAnkle,nHip]);
GA.Sim.Con.MinSat = [-maxAnkle,-maxHip];
GA.Sim.Con.MaxSat = [ maxAnkle, maxHip];

% Simulation parameters
GA.Sim.IC = [start_slope, start_slope, 0, 0, zeros(1, GA.Sim.Con.stDim)];
GA.Sim = GA.Sim.SetTime(0,0.03,20);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;

GA.FitFcn = {1, @MOOGA.VelFit;
             2, @MOOGA.NrgEffFit;
             3:10, @MOOGA.VelRangeFit;
             11, @MOOGA.EigenFit};
GA.FitIDs = [1,2,3]; % Velocity and average COT
GA.FitMinMax = [1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1];

% GA.FitFcn = {1, @MOOGA.VelFit;
%              2, @MOOGA.EigenFit;
%              3:10, @MOOGA.VelRangeFit};
% GA.FitIDs = [3,4]; % Velocity range and average COT
% GA.FitMinMax = [1, 1, 1, 1, -1, 1, -1, 1, -1, 1];
% % Pareto always looks to maximize the fitness value but sometimes we want
% % to check values that don't go into the pareto selection and we want to
% % minimize (or get the maximum negative value, e.g. a negative slope)

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

end
