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
Sim.Con = Sim.Con.HandleExtFB(Sim.IC(Sim.ModCo),Sim.IC(Sim.ConCo));

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
for f = 1:GA.NFit
    thisFit(f) = GA.FitFcn{f}(Sim);
end
disp(['Genome ',num2str(ID),' results: ',num2str(thisFit)]);

end

