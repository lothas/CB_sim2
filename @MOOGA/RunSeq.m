function RunSeq( GA, varargin )
%RUNSEQ Runs a simulation with the selected genome

tend = 45;
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

% GA.FitFcn{end+1} = @MOOGA.ZMPFit;
% GA.NFit = length(GA.FitFcn);

Sim = deepcopy(GA.Sim);
Sim.Graphics = 1;
Sim.EndCond = 2; % Run until converge

Sim = GA.Gen.Decode(Sim, GA.Seqs(ID,:,Generation));
disp(GA.Gen.seq2str(GA.Seqs(ID,:,Generation)));
% Simulation parameters
Sim = Sim.SetTime(0,0.05,tend);

% Set internal parameters (state dimensions, events, etc)
Sim = Sim.Init();

% Some more simulation initialization
Sim.Mod.LegShift = Sim.Mod.Clearance;
Sim.Con = Sim.Con.Reset(Sim.IC(Sim.ConCo));
Sim.Con.FBType = 2;
% Sim.Con = Sim.Con.HandleExtFB(Sim.IC(Sim.ModCo),...
%                 Sim.IC(Sim.ConCo),Sim.Env.SurfSlope(Sim.Mod.xS));

% Simulate
Sim = Sim.Run();

% convT = Sim.Out.T(end);
% return

% Calculate the genome's fitness
thisFit = zeros(1,max(cell2mat(GA.FitFcn(:,1)')));
thisOuts = cell(1,GA.NFit);
for f = 1:GA.NFit
    % Preprocessing for ZMPFit
    if ~isempty(strfind(func2str(GA.FitFcn{f,2}),'ZMPFit'))
        % Prepare all the required vectors
        % (torques, state, etc) and put them in wSim.Out

        % ZMP Fit should be the last one as it uses the
        % output from all previous fitness functions
        Outs = find(~cellfun(@isempty,thisOuts));
        for o = 1:length(Outs)
            Sim.Out = Sim.JoinOuts(thisOuts{Outs(o)});
        end
        Sim.Con.FBType = 2;
    end
    
    FitInd = GA.FitFcn{f,1};
    try
        [thisFit(FitInd),thisOuts{f}] = GA.FitFcn{f,2}(Sim);
    catch
        FuncName = MOOGA.GetFitFcnName(GA.FitFcn{f,2});
        switch FuncName
            case 'VelFit'
                [thisFit(FitInd),thisOuts{f}] = MOOGA.VelFit(Sim);
            case 'NrgEffFit'
                [thisFit(FitInd),thisOuts{f}] = MOOGA.NrgEffFit(Sim);
            case 'EigenFit'
                [thisFit(FitInd),thisOuts{f}] = MOOGA.EigenFit(Sim);
            case 'UphillFitRun'
                [thisFit(FitInd),thisOuts{f}] = MOOGA.UphillFitRun(Sim);
            case 'DownhillFitRun'
                [thisFit(FitInd),thisOuts{f}] = MOOGA.DownhillFitRun(Sim);
            case 'ZMPFit'
                [thisFit(FitInd),thisOuts{f}] = MOOGA.ZMPFit(Sim);
            otherwise
                thisFit(FitInd) = 0;
                thisOuts{f} = [];
        end
    end
end
disp(['Genome ',num2str(ID),' results: ',num2str(thisFit)]);

end

