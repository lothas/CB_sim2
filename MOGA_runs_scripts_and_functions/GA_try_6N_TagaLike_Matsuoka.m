function [  ] =...
    GA_try_6N_TagaLike_Matsuoka(whichCase,fileIn,trainingDataFile)
% Run MOOGA using only Vel, Nrg and Speed-range fitness but limiting the
% simulation with a bounded foot size (ZMP threshold)


GA = MOOGA(20,500);
GA = GA.SetFittest(15,15,0.5);
GA.JOAT = 2; GA.Quant = 0.7;

GA.FileIn = fileIn;

FileName_start = 'VGAM_6N_TagaLike_';
FileName_date = datestr(now,'mm_dd_hh_MM');
FileName_extra = '_same_tonicInputs_20Gen_500Genes_';

% Check ankle Torque?
        %   Assign a check function to check if a CPG doesn't have an oscillatory
        %   Ankle joint.
GA.genomeChevkFcn = [];%@genomeChevkFcn;

switch whichCase
    case 'GA only'
        % Use NN?
        use_NN = 0;
%         out file name:
        GA.FileOut = [FileName_start,FileName_date,FileName_extra,...
            '_GA_only','.mat'];

   case 'GA + rescale'
        use_NN = 0;
        GA.rescaleFcn = @rescaleFcn;
        GA.FileOut = [FileName_start,FileName_date,FileName_extra,...
            '_rescale_only','.mat']; 
        
    case 'GA + NN_reg'
        use_NN = 'NN_reg';
        GA.FileOut = [FileName_start,FileName_date,FileName_extra,...
            '_NN_reg_only','.mat'];
        
    case 'GA + NN_classi'
        use_NN = 'NN_classi';
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
load('MatsuokaGenome_4Neuron_tagaLike.mat','Keys','Range','N',...
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

        inFilenames = {trainingDataFile};
        
        MML.sample_genes = {'\tau_r','beta','4neuron_taga_like'}; 
        MML.target_genes = {'period'};

        [samples, targets, normParams] = MML.prepare_reg_NNData('6N_CPG',inFilenames, maxN);
        MML.normParams = normParams;

        architecture = [10];
        [net, tr] = ...
                MML.train_reg_NN(samples, targets, architecture, NNSamples);
        save(GANN_file,'net','tr');

        GA.NN_reg = net;
        GA.NN_reg_Fcn = @NN_reg_Fcn;
    case {'NN_classification','NN_classi'}
        GANN_file = 'MatsGANN_classi.mat';

        maxN = 250000;
        
        inFilenames = {trainingDataFile};%
%         MML.sample_genes = {'\tau_r','beta','6neuron_taga_like'}; 
        MML.sample_genes = {'\tau_r','beta','6neuron_taga_like_symm'}; 
        MML.target_genes = {'n_osc and osc classes'};
        [samples, targets] = ...
            MML.prepare_classi_NNData('6N_CPG',inFilenames, maxN);

        MML.normParams = [];

        architecture = [10];
        [net, ~] = MML.train_classi_NN(samples, targets, architecture);
        save(GANN_file,'net');

        GA.NN_classi = net;
        GA.NN_classi_Fcn = @NN_classi_Fcn;
end

function seq = NN_reg_Fcn(Gen, net, seq, X, T)

    % % don't need to run Sim again. run the NN on all CPG's??
    
    [~, periods, ~, ~, ~] = MML.processResults(X, T);
    % don't do anything if CPG IS oscillating
    if ~any(isnan(periods))
        return
    end

    % Use NN to select best value for tau gene
    desPeriod = MML.perLim(1) + ...
                 rand()*(MML.perLim(2)-MML.perLim(1));

    seq = MML.get_reg_NNPar(net, seq, desPeriod);
    
    ids = seq < Gen.Range(1,:) | seq > Gen.Range(2,:);
    seq(ids) = min(max(seq(ids),Gen.Range(1,ids)),Gen.Range(2,ids));
end

%%
function noGO_flag = genomeChevkFcn(X,T)
    % If the ankle torque is a posivite constant (larger than zero and with
    % no period) than don't use this genome.
    [~, periods, signals, ~, ~] = MML.processResults(X, T);
    
    % check that the ankle torque is positive:
    ind2analize = ceil(0.7*size(signals,2));
    lastOfAnkleSignal = signals(1,ind2analize:end);
    
    %  dont take CPGs if the ankle torque is not oscillating and positive
    if (mean(lastOfAnkleSignal) > 0.001) && isnan(periods(1,1))
        noGO_flag = 1;
    else % if the hip torque is constant than ignore it as well
        if isnan(periods(2,1))
            noGO_flag = 1;
        else
            noGO_flag = 0;
        end
    end
    
end
%% % % % First Class_NN vesion:
% % Ruffle a sample from a random group.
% 
% function seq = NN_classi_Fcn(Gen, net, seq, X, T)
%     
%     [~, periods, ~, ~, ~] = MML.processResults(X, T);
%     % don't do anything if CPG IS oscillating
%     if ~any(isnan(periods)) 
%         return
%     end
% 
%     rand_seq = MML.Gen.RandSeq(5000);
%     seq = MML.get_classi_NNPar(net, rand_seq);
% 
%     ids = seq < Gen.Range(1,:) | seq > Gen.Range(2,:);
%     seq(ids) = min(max(seq(ids),Gen.Range(1,ids)),Gen.Range(2,ids));
% end
%% % Improve #1 Class_NN version:
% % Ruffle it from the cross mutated version of the last generation
function seq = NN_classi_Fcn(Gen,net,seq,lastGen,lastGenes, X, T)
    % get MOGA class instead of Gen class. Instead off ruffling random 
    %   CPGs, it's getting a random CPGs from the Mutated versions of the
    %   topPop.
    
    [~, periods, ~, ~, ~] = MML.processResults(X, T);
    % don't do anything if CPG IS oscillating
    if ~any(isnan(periods)) 
        return
    end
    
    % Give 10% chance for a random good sample to get selected (incresing
    % exploration):
    if randi(10) == 10
        rand_seq = MML.Gen.RandSeq(100000);
        seq = MML.get_classi_NNPar(net, rand_seq);
    else
        try % % try to select a mutated version of the topPop:
            rand_seq = [lastGenes;...
                MML.Gen.RandSeq(2000)];
            seq = MML.get_classi_NNPar(net, rand_seq);
        catch % if there are no good CPGs ruffle moere random samples
            rand_seq = MML.Gen.RandSeq(100000);
            seq = MML.get_classi_NNPar(net, rand_seq);
        end
    end

    ids = seq < Gen.Range(1,:) | seq > Gen.Range(2,:);
%     seq(ids) = min(max(seq(ids),Gen.Range(1,ids)),Gen.Range(2,ids));
end

%% % Improve #2 Class_NN version:
% % % Ruffle it from the cross mutated version of the last generation
% % % ALSO, anable 10% of getting not making any change at all.
% function seq = NN_classi_Fcn(Gen,net,seq,lastGen,lastGenes, X, T)
%     % get MOGA class instead of Gen class. Instead off ruffling random 
%     %   CPGs, it's getting a random CPGs from the Mutated versions of the
%     %   topPop.
%     
%     % get radom integer from '1' to '10' if this number is 10 than dont use
%     % the NN.
%     if randi(10) > 9
%         return;
%     end
%     
%     [~, periods, ~, ~, ~] = MML.processResults(X, T);
%     % don't do anything if CPG IS oscillating
%     if ~any(isnan(periods)) 
%         return
%     end
%     
%     try
%         rand_seq = [lastGenes;...
%             MML.Gen.RandSeq(2000)];
%         seq = MML.get_classi_NNPar(net, rand_seq);
%     catch % if there are no good CPGs
%         rand_seq = MML.Gen.RandSeq(100000);
%         seq = MML.get_classi_NNPar(net, rand_seq);
%     end
% 
%     ids = seq < Gen.Range(1,:) | seq > Gen.Range(2,:);
%     % seq(ids) = min(max(seq(ids),Gen.Range(1,ids)),Gen.Range(2,ids));
% end
%% % Rescaling function:
function seq = rescaleFcn(Gen, seq, X, T)
    [~, periods, ~, ~, ~] = MML.processResults(X, T);

    % don't do anything if CPG is not oscillating
    if ~any(isnan(periods))
        return
    end
    inputPeriod = mean(periods);

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
GA.Sim.Con = GA.Sim.Con.SetAnkles_selection(2); % set the simulation to work with two ankle joints
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

% GA.FitFcn = {1, @MOOGA.VelFit};
% GA.FitIDs = [1]; % Velocity and average COT
% GA.FitMinMax = [1];

% GA.FitFcn = {1, @MOOGA.VelFit;
%              2, @MOOGA.NrgEffFit;
%              3, @MOOGA.EigenFit};
% GA.FitIDs = [1,2,3]; % Velocity and average COT
% GA.FitMinMax = [1, 1, 1];

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