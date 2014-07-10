function RunSeq( GA, varargin )
%RUNSEQ Runs a simulation with the selected genome

tend = 'inf';
switch nargin 
    case 2
        Generation = GA.Progress;
        ID = varargin{1};
    case 3
        Generation = varargin{1};
        ID = varargin{2};
    case 4
        Generation = varargin{1};
        ID = varargin{2};
        tend = varargin{3};
    otherwise
        Generation = GA.Progress;
        TopIDs = GA.GetTopPop(GA.Fittest(1));
        ID = randsample(TopIDs,1);
end

Sim = deepcopy(GA.Sim);
Sim.Graphics = 1;
Sim.EndCond = 2; % Run until converge

Sim = GA.Gen.Decode(Sim, GA.Seqs(ID,:,Generation));
disp(num2str(GA.Seqs(ID,:,Generation)',5));
% Simulation parameters
Sim = Sim.SetTime(0,0.09,tend);

% Set internal parameters (state dimensions, events, etc)
Sim = Sim.Init();

% Some more simulation initialization
Sim.Mod.LegShift = Sim.Mod.Clearance;
% Sim.Con = Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));
Sim.Con.FBType = 0;
% Sim.Con = Sim.Con.HandleExtFB(Sim.IC(Sim.ModCo),Sim.IC(Sim.ConCo));

% Simulate
Sim = Sim.Run();

% Calculate eigenvalues
if Sim.Out.Type == 5
    [EigVal,EigVec] = Sim.Poincare();
    % Do some plots
    disp(EigVal);
else
    EigVal = 2*ones(4,1);
    disp(Sim.Out.Text);
end

% Calculate the genome's fitness
thisFit = zeros(1,GA.NFit);
thisOuts = cell(1,GA.NFit);
for f = 1:GA.NFit
    % Preprocessing for ZMPFit
    if strcmp(func2str(GA.FitFcn{f}),...
        '@(varargin)GA.ZMPFit(varargin{:})')
        % Prepare all the required vectors
        % (torques, state, etc) and put them in wSim.Out

        % ZMP Fit should be the last one as it uses the
        % output from all previous fitness functions
        Outs = find(~cellfun(@isempty,thisOuts));
        X = []; T = 0;
        SuppPos = []; Torques = [];
        for o = 1:length(Outs)
            X = [X;thisOuts{Outs(o)}.X];
            T = [T;thisOuts{Outs(o)}.T+T(end)];
            SuppPos = [SuppPos;thisOuts{Outs(o)}.SuppPos];
            Torques = [Torques;thisOuts{Outs(o)}.Torques];
        end
        Sim.Out.X = X;
        Sim.Out.T = T(2:end);
        Sim.Out.SuppPos = SuppPos;
        Sim.Out.Torques = Torques;
        Sim.Con.FBType = 2;
    end
    
    [thisFit(f),thisOuts{f}] = GA.FitFcn{f}(Sim);
end
disp(['Genome ',num2str(ID),' results: ',num2str(thisFit)]);

end

