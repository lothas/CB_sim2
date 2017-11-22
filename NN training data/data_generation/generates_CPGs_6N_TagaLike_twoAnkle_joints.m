
close all; clc; clear all;

%% Run a more general case of Taga-like CPG (six types of weights)
generate_GenomeFile('6N_tagaLike_2Ank_torques')

% Set up the genome
load('MatsuokaGenome_4Neuron_tagaLike.mat','N');

% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 2*N;
clear N

header = sprintf('tau ratio is equal to %d \n',MML.Sim.Con.tau_ratio);
header = [header,sprintf('data is for 6N TagaLike case with same C for each joint\n')];
header = [header,sprintf('seq Order: \n')];
header = [header,sprintf('"tau","b","c","w1","w2","w3","w4","w5","w6","k_tau","k_c" \n')];
header = [header,sprintf('"tau" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,1),MML.Gen.Range(2,1))];
header = [header,sprintf('"b" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,2),MML.Gen.Range(2,2))];
header = [header,sprintf('"c" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,3),MML.Gen.Range(2,3))];
header = [header,sprintf('"w_1" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,4),MML.Gen.Range(2,4))];
header = [header,sprintf('"w_2" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,5),MML.Gen.Range(2,5))];
header = [header,sprintf('"w_3" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,6),MML.Gen.Range(2,6))];
header = [header,sprintf('"w_4" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,7),MML.Gen.Range(2,7))];
header = [header,sprintf('"w_5" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,8),MML.Gen.Range(2,8))];
header = [header,sprintf('"w_6" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,9),MML.Gen.Range(2,9))];
header = [header,sprintf('"k_tau" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,10),MML.Gen.Range(2,10))];
header = [header,sprintf('"k_c" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,11),MML.Gen.Range(2,11))];

disp(header);

%% %% Run a more specific case of Taga-like CPG (four types of weights)
% in this case there is more symmetry in the CPG
generate_GenomeFile('6N_tagaLike_2Ank_torques_symm');

% Set up the genome
load('MatsuokaGenome_4Neuron_tagaLike.mat','N');

% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;
MML.nNeurons = 2*N;
clear N

header = sprintf('tau ratio is equal to %d \n',MML.Sim.Con.tau_ratio);
header = [header,sprintf('data is for 6N TagaLike case with same C for each joint\n')];
header = [header,sprintf('seq Order: \n')];
header = [header,sprintf('"tau","b","c","w1","w2","w3","w4","k_tau","k_c" \n')];
header = [header,sprintf('"tau" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,1),MML.Gen.Range(2,1))];
header = [header,sprintf('"b" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,2),MML.Gen.Range(2,2))];
header = [header,sprintf('"c" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,3),MML.Gen.Range(2,3))];
header = [header,sprintf('"w_1" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,4),MML.Gen.Range(2,4))];
header = [header,sprintf('"w_2" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,5),MML.Gen.Range(2,5))];
header = [header,sprintf('"w_3" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,6),MML.Gen.Range(2,6))];
header = [header,sprintf('"w_4" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,7),MML.Gen.Range(2,7))];
header = [header,sprintf('"k_tau" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,8),MML.Gen.Range(2,8))];
header = [header,sprintf('"k_c" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,9),MML.Gen.Range(2,9))];

disp(header);
%% Train data: get All osc and n-osc (train data for classifier)
N = 1000; % the number of samples
CPGs_num = 0;
round_count = 0;
max_round = 1000; % maximum iteraion for while loop (saftey reasons:)
results = [];

wanted_num_CPGs = 200000;

disp('start with the sim:');

while (CPGs_num < wanted_num_CPGs) && (round_count <= max_round)
    
    rand_seq = MML.Gen.RandSeq(N);
    parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
%     for i=1:N
%         disp(['at sim #',num2str(i)]);
        [out, ~, signal] = MML.runSim(rand_seq(i,:));
            % Prepare output:
        % Parameters
        results_temp(i).seq = rand_seq(i,:);

        % Results- caculate perdiods using different methods:
        results_temp(i).periods = out.periods;
        
        results_temp(i).x0 = out.x0;
        results_temp(i).neuronActive = out.neuronActive;
        results_temp(i).neuronOsc = out.neuronOsc;
    end 

    results = [results,results_temp];
    CPGs_num = length(results);

    disp(['so far we have ',num2str(CPGs_num),' CPGs'])

    round_count = round_count+1;
    
    clear results_temp rand_seq

end

disp('sim end...');

save('MatsRandomRes_6N_TagaLike_TrainingSet_2.mat',...
    'results','header');

clear N

%% Train data: get Only osc ( for regression NN)
N = 1000; % the number of samples
CPGs_num = 0;
round_count = 0;
max_round = 1000; % maximum iteraion for while loop (saftey reasons:)
results = [];

wanted_num_CPGs = 20000;

disp('start with the sim:');

while (CPGs_num < wanted_num_CPGs) && (round_count <= max_round)
    
    rand_seq = MML.Gen.RandSeq(N);
    parfor i=1:N % Simulate and calculate the frequecy (also calc from Matsuoka extimation)
%     for i=1:N
%         disp(['at sim #',num2str(i)]);
        [out, ~, signal] = MML.runSim(rand_seq(i,:));
            % Prepare output:
        % Parameters
        results_temp(i).seq = rand_seq(i,:);

        % Results- caculate perdiods using different methods:
        results_temp(i).periods = out.periods;
        
        results_temp(i).x0 = out.x0;
        results_temp(i).neuronActive = out.neuronActive;
        results_temp(i).neuronOsc = out.neuronOsc;
    end 
    
    periods = horzcat(results_temp.periods);
    
    % Filter CPG's where not both signals oscillating:
    osc_ids = ~isnan(periods);
    osc_ids = osc_ids(1,:) & osc_ids(2,:) & osc_ids(3,:);

    % Filter CPG's where the is a big difference between hip and ankle:
    periods_ratios = (periods(1,:)./periods(2,:));
    diff_ids1 = (periods_ratios >  0.85) & (periods_ratios <  1.15); 

    periods_ratios = (periods(1,:)./periods(3,:));
    diff_ids2 = (periods_ratios >  0.85) & (periods_ratios <  1.15); 

    good_ids = osc_ids & diff_ids1 & diff_ids2;
    
    results = [results,results_temp(good_ids)];
    CPGs_num = length(results);

    disp(['so far we have ',num2str(CPGs_num),' CPGs'])

    round_count = round_count+1;
    
    clear results_temp rand_seq

end

disp('sim end...');

save('MatsRandomRes_6N_TagaLike_TrainingSet_RegressionNN.mat',...
    'results','header');

clear N
%% Run CB:

% Set up the genome
load('MatsuokaGenome_4Neuron_tagaLike.mat','Keys','Range','N',...
    'nAnkle','nHip','maxAnkle', 'maxHip','MutDelta0','MutDelta1');

GA = MOOGA(2,100);
GA.Gen = Genome(Keys, Range);

GA.Sim = Simulation();
GA.Sim.Graphics = 1;
GA.Sim.EndCond = 2; % Run until converge (or fall)

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

GA.NFit = size(GA.FitIDs,2);
GA.Sim.PMFull = 1; % Run poincare map on all coords

GA = GA.InitGen();

thisSeq = results(n).seq;

GA.RunSeq(thisSeq)


clear Keys Range N nAnkle nHip maxAnkle...
    maxHip MutDelta0 MutDelta1 GA Sim
