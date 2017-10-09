function [  ] = GA_try_2N_Symm_Matsuoka(whichCase,fileIn)
% Run MOOGA using only Vel, Nrg and Speed-range fitness but limiting the
% simulation with a bounded foot size (ZMP threshold)



GA = MOOGA(20,500);
GA = GA.SetFittest(15,15,0.5);
GA.JOAT = 2; GA.Quant = 0.7;

GA.FileIn = fileIn;

FileName_start = 'VGAM_2N_symm_';
FileName_date = datestr(now,'mm_dd_hh_MM');
FileName_extra = '_larger_range';

switch whichCase
    case 'GA only'
        % Use NN?
        use_NN = 0;
        % Rescale?
        GA.rescaleFcn = [];
%         out file name:
        GA.FileOut = ['VGAM_2N_symm_',datestr(now,'mm_dd_hh_MM'),...
            '_GA_only','.mat'];

   case 'GA + rescale'
        use_NN = 0;
        GA.rescaleFcn = @rescaleFcn;
        GA.FileOut = [FileName_start,FileName_date,FileName_extra,...
            '_rescale_only','.mat']; 
        
    case 'GA + NN_reg'
        use_NN = 'NN_reg';
        GA.rescaleFcn = [];
        GA.FileOut = [FileName_start,FileName_date,FileName_extra,...
            '_NN_reg_only','.mat'];
        
    case 'GA + NN_classi'
        use_NN = 'NN_classi';
        GA.rescaleFcn = [];
        GA.FileOut = [FileName_start,FileName_date,FileName_extra,...
            '_NN_classi_only','.mat'];
        
    
    case 'GA + NN_reg + rescale'
        use_NN = 'NN_reg';
        GA.rescaleFcn = @rescaleFcn;
        GA.FileOut = [FileName_start,FileName_date,FileName_extra,...
            '_NNreg_and_rescale','.mat'];
        
    case 'GA + NN_classi + rescale'
        use_NN = 'NN_classi';
        GA.rescaleFcn = @rescaleFcn;
        GA.FileOut = [FileName_start,FileName_date,FileName_extra,...
            '_NNclassi_and_rescale','.mat'];
end
        

GA.Graphics = 0;
GA.ReDo = 1;

% Set up the genome
load('MatsuokaGenome_2Neuron_Symm.mat','Keys','Range','N',...
    'nAnkle','nHip','maxAnkle', 'maxHip','MutDelta0','MutDelta1');

GA.Gen = Genome(Keys, Range);

MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.nNeurons = 2*N;

% Use NN?
switch use_NN
    case {'NN_regression','NN_reg'}
        GANN_file = 'MatsGANN_reg.mat';

        maxN = 250000;
        NNSamples = 500;

        inFilenames =...
            {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_1.mat'};
        
        MML.sample_genes = {'\tau_r','2neuron_symm_weights'}; 
        MML.target_genes = {'beta'};

    %     MML.sample_genes = {'2neuron_symm_weights'}; 
    %     MML.target_genes = {'\tau_r','beta'};

        [samples, targets, normParams] = MML.prepare_reg_NNData('2N_CPG',inFilenames, maxN);
        MML.normParams = normParams;

    %     if exist(GANN_file, 'file') ~= 2
            architecture = [10];
            [net, tr] = ...
                    MML.train_reg_NN(samples, targets, architecture, NNSamples);
            save(GANN_file,'net','tr');
    %     else
    %         GANN_net = load(GANN_file);
    %         net = GANN_net.net;
    %     end

        GA.NN_reg = net;
        GA.NN_reg_Fcn = @NN_reg_Fcn;
    case {'NN_classification','NN_classi'}
        GANN_file = 'MatsGANN_classi.mat';

        maxN = 250000;

%         inFilenames =...
%             {'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_All_1.mat',...
%             'MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_All_2.mat'};

        inFilenames = {'MatsRandomRes_2Neurons_symm_very_large_b_and_tau_and_a_All_1.mat'};
        
        MML.sample_genes = {'\tau_r','beta','2neuron_symm_weights'}; 
        MML.target_genes = {'n_osc and osc classes'};

        [samples, targets, normParams] = ...
            MML.prepare_classi_NNData('2N_CPG',inFilenames, maxN);
        MML.normParams = normParams;

    %     if exist(GANN_file, 'file') ~= 2
            architecture = [10];
            [net, ~] = MML.train_classi_NN(samples, targets, architecture);
            save(GANN_file,'net');
    %     else
    %         GANN_net = load(GANN_file);
    %         net = GANN_net.net;
    %     end

        GA.NN_classi = net;
        GA.NN_classi_Fcn = @NN_classi_Fcn;
end

function seq = NN_reg_Fcn(Gen, net, seq, X, T)

    % % don't need to run Sim again. run the NN on all CPG's??
    
    [~, periods, ~, ~, ~] = MML.processResults(X, T);
    % don't do anything if CPG IS oscillating
    % in the case with 2N (actuated hip only) check the hip_period only:
    if ~isnan(periods(2)) 
        return
    end

    % Use NN to select best value for tau gene
    desPeriod = MML.perLim(1) + ...
                 rand()*(MML.perLim(2)-MML.perLim(1));

    seq = MML.get_reg_NNPar(net, seq, desPeriod);
    
    ids = seq < Gen.Range(1,:) | seq > Gen.Range(2,:);
    seq(ids) = min(max(seq(ids),Gen.Range(1,ids)),Gen.Range(2,ids));
end

function seq = NN_classi_Fcn(Gen, net, seq, X, T)
    
    [~, periods, ~, ~, ~] = MML.processResults(X, T);
    % don't do anything if CPG IS oscillating
    % in the case with 2N (actuated hip only) check the hip_period only:
    if ~isnan(periods(2)) 
        return
    end

    rand_seq = MML.Gen.RandSeq(500);
    seq = MML.get_classi_NNPar(net, rand_seq);

    ids = seq < Gen.Range(1,:) | seq > Gen.Range(2,:);
    seq(ids) = min(max(seq(ids),Gen.Range(1,ids)),Gen.Range(2,ids));
end

function seq = rescaleFcn(Gen, seq, X, T)
    [~, periods, ~, ~, ~] = MML.processResults(X, T);

    % don't do anything if CPG is not oscillating
    % in the case with 2N (actuated hip only) check the hip_period only:
    if isnan(periods(2)) 
        return
    end
    inputPeriod = periods(2);

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
