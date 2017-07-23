function [avg_sim_time,simOutType] = run_CPG_with_CB(results)
% run Matsuoka CPG's with the CB model to check simulation time

% Set up the genome
genome_file = 'MatsuokaGenome.mat';
load(genome_file);

MML = MatsuokaML();
MML.perLim = [0.68 0.78];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.nNeurons = 2*N;

% Set up the simulations
Sim = Simulation();
Sim.Graphics = 0;
% Sim.EndCond = 2; % Run until converge (or fall)
Sim.EndCond = 0; % Run until the end of time

% Set up the compass biped model
Sim.Mod = Sim.Mod.Set('I',0,'damp',0,'A2T',0.16,'A2H',0.12);

% Set up the terrain
start_slope = 0;
Sim.Env = Sim.Env.Set('Type','inc','start_slope',start_slope);

% Initialize the controller
Sim.Con = Matsuoka;
Sim.Con.startup_t = 1.0; % Give some time for the neurons to converge
                            % before applying a torque
Sim.Con.FBType = 0; % no slope feedback
Sim.Con.nPulses = N;
Sim.Con.stDim = 4*N;
Sim.Con = Sim.Con.SetOutMatrix([nAnkle,nHip]);
Sim.Con.MinSat = [-maxAnkle,-maxHip];
Sim.Con.MaxSat = [ maxAnkle, maxHip];

% Simulation parameters
Sim.IC = [start_slope, start_slope, 0, 0, zeros(1, Sim.Con.stDim)];
Sim = Sim.SetTime(0,0.05,15);

% Some more simulation initialization
Sim.Mod.LegShift = Sim.Mod.Clearance;

seq = vertcat(results(:).seq);
tic
n = length(results);

t_start = tic;
conv_counter = 0;
t_conv_sim = 0;
simOutType = zeros(1,n);

for i=1:n
    t_cur = tic;
    trySim = deepcopy(Sim);
    trySim = MML.Gen.Decode(Sim, seq(i,:));
    trySim.Con.s_in = 0;
    trySim.Con = trySim.Con.Adaptation();
    trySim = trySim.SetTime(0,0.01,10);
    trySim = trySim.Init();
    trySim.Run();
    
    simOutType(1,i) = trySim.Out.Type;
    
    % check if the simulation convarged:
    if trySim.Out.Type==6 % any(trySim.Out.Type == [0,5,6])
        % Out.Type==0 : reached the end of tSpan
        t_conv_sim = t_conv_sim + toc(t_cur);
        conv_counter = conv_counter + 1;
    end
    
    disp(['at i = ',num2str(i),...
        ' conv_count = ',num2str(conv_counter),...
        ' ourType = ',num2str(trySim.Out.Type)]);
end
t_elapsed = toc(t_start);
avg_sim_time = t_elapsed/n;
disp(['avg sim time for all simulations is ',num2str(avg_sim_time),' [sec]']);

disp(['avg sim time for conv simulations is ',num2str(t_conv_sim/conv_counter),' [sec]']);
end
