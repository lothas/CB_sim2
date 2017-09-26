
close all; clc; clear all;

%%
genome_file = 'MatsuokaGenome_2Neuron_Symm.mat';
nAnkle = 1;%1; % Number of ankle torques
nHip = 1;   % Number of hip torques
maxAnkle = 10;   % Max ankle torque
maxHip = 10;    % Max hip torque
Mamp = [10,10];
mamp = 0*Mamp;
N = nAnkle+nHip;
% %     % 2neuron symmetric specific range%%

       % % % Narrow b Narrow W Narrow tau
Mw = 5*1;
mw = 0*Mw;
Keys = {'\tau_r', 'beta','amp_2n_same_inputs',    '2neuron_symm_weights', 'ks_\tau',     'ks_c_2n_symm', 'IC_matsuoka';
              1 ,      1,                   2,                         1,        1 ,          1,            0 };
Range = {  0.02 ,    0.2,               [0,0],                         0,   -0.001 ,       -0.2; % Min
           0.10  ,   2.5,             [10,10],                         5,    0.001 ,       0.2}; % Max

       
MutDelta0 = 0.04;   MutDelta1 = 0.02;

save(genome_file, 'nAnkle', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

clear all
%%
% define the class:
MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.tStep = 0.05;
MML.tEnd = 15;

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


header = sprintf('tau ratio is equal to %d \n',MML.Sim.Con.tau_ratio);
header = [header,sprintf('data is for 2N symmetric case \n')];
header = [header,sprintf('seq Order: \n')];
header = [header,sprintf('"tau","b","c","NR","a" \n')];
header = [header,sprintf('"tau" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,1),MML.Gen.Range(2,1))];
header = [header,sprintf('"b" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,2),MML.Gen.Range(2,2))];
header = [header,sprintf('"c" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,3),MML.Gen.Range(2,3))];
header = [header,sprintf('"a" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,5),MML.Gen.Range(2,5))];


save('MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_All_1.mat',...
    'results','header');

clear N

%% Train data: (get good ones)
N = 1000; % the number of samples
good_CPGs_num = 0;
round_count = 0;
max_round = 1000; % maximum iteraion for while loop (saftey reasons:)
results = [];

wanted_num_CPGs = 1000;

disp('start with the sim:');

while (good_CPGs_num < wanted_num_CPGs) && (round_count <= max_round)
    rand_seq = MML.Gen.RandSeq(N);
%     results_temp(N) = [];
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

    % get the periods
    periods = horzcat(results_temp(:).periods);
    periods = periods(2,:);
    
    % get oscillatory CPGs ids:
    osc_ids = ~isnan(periods);
    
%     % get oscillatory CPGs ids in period range::
%     osc_ids = ~isnan(periods) &...
%         ((periods > MML.perLimOut(1,1)) & (periods < MML.perLimOut(1,2)));
    
    % get neural oscillation check:
    neuronOsc = (vertcat(results_temp(:).neuronOsc))';
    neuronOsc = neuronOsc(3:4,:);   % only take hip neurons
    osc_check_ids = true(1,N);
    for k=1:size(neuronOsc,1) % run on all neurons
        % mark as "good" if we have at least one osc neuron
        osc_check_ids = osc_check_ids & ~neuronOsc(k,:);
    end
    
    good_ids = osc_ids & ~osc_check_ids;

    str1 = sprintf('round %d :    osc_CPGs: %d , ',round_count,sum(osc_ids));
    str1 = [str1,sprintf('    osc_MN_ids: %d ,    ',sum(osc_check_ids))];
    str1 = [str1,sprintf('    tot_CPG: %d \n',sum(good_ids))];    
    disp(str1);
    
    results = [results,results_temp(good_ids)];
    good_CPGs_num = length(results);

    disp(['so far we have ',num2str(good_CPGs_num),' CPGs'])

    round_count = round_count+1;
    
    clear results_temp periods rand_seq osc_check_ids neuronOsc
    clear good_ids osc_ids str1

end

disp('sim end...');


header = sprintf('tau ratio is equal to %d \n',MML.Sim.Con.tau_ratio);
header = [header,sprintf('data is for 2N symmetric case \n')];
header = [header,sprintf('seq Order: \n')];
header = [header,sprintf('"tau","b","c","NR","a" \n')];
header = [header,sprintf('"tau" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,1),MML.Gen.Range(2,1))];
header = [header,sprintf('"b" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,2),MML.Gen.Range(2,2))];
header = [header,sprintf('"c" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,3),MML.Gen.Range(2,3))];
header = [header,sprintf('"a" in range ( %.2f , %.2f ) \n',...
    MML.Gen.Range(1,5),MML.Gen.Range(2,5))];


save('MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_osc_inRange_test_group1.mat',...
    'results','header');
% save('MatsRandomRes_2Neurons_symm_Narrow_b_Narrow_W_Narrow_tau_only_osc_2.mat',...
%     'results','header');

clear N
%% plot CPG output:
close all;

% n=18;
n = randsample(1:length(results),1);

[out, ~, signal] = MML.runSim(results(n).seq);
figure;
subplot(2,1,1);
plot(signal.T,signal.X);
xlabel('time[sec]');    ylabel('X_i');
legend('x1','x2','x3','x4','xp1','xp2','xp3','xp4');
title({'X_i over time',...
    ['id #',num2str(n),...
    '    periods: ',...
    num2str(out.periods')],...
    ['PerOK1: ',num2str(out.perOK1),...
    '    PerOK2: [',num2str(out.perOK2'),']'],...
    ['neuronActive: ',num2str(out.neuronActive),...
    '    neuronOsc: [',num2str(out.neuronOsc),']']});
subplot(2,1,2)
plot(signal.T,signal.signal(1,:),'b',signal.T,signal.signal(2,:),'r');
legend('ankle','hip');
clear signal out


%% Run CB:

% Set up the genome
load('MatsuokaGenome_2Neuron_symm.mat','Keys','Range','N',...
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

thisSeq = GA_temp.Seqs(1,:,5);
% thisSeq = results(n).seq;

GA.RunSeq(thisSeq)


clear Keys Range N nAnkle nHip maxAnkle...
    maxHip MutDelta0 MutDelta1 GA Sim
