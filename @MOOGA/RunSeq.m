function varargout = RunSeq( GA, varargin )
%RUNSEQ Runs a simulation with the selected genome

switch nargin 
    case 2
        if length(varargin{1}) == 1
            Generation = GA.Progress;
            ID = varargin{1};
        else
            thisSeq = varargin{1};
        end
    case 3
        Generation = varargin{1};
        ID = varargin{2};
    case 4
        Generation = varargin{1};
        ID = varargin{2};
        tend = varargin{3};
        GA.Sim = GA.Sim.SetTime(GA.Sim.tstart, ...
                                GA.Sim.tstep, tend);
    otherwise
        Generation = GA.Progress;
        TopIDs = GA.GetTopPop(GA.Fittest(1));
        ID = randsample(TopIDs,1);
end

if ~exist('thisSeq','var')
    thisSeq = GA.Seqs(ID,:,Generation);
end

% Set-up the simulation
wSim = deepcopy(GA.Sim);
    
wSim = GA.Gen.Decode(wSim,thisSeq);

if nargout == 0
    wSim.Graphics = 1;
    wSim.EndCond = 2; % Run until converge
    disp(GA.Gen.seq2str(thisSeq));
end

% Set internal parameters (state dimensions, events, etc)
wSim = wSim.Init();

% start_slope = -0.2*pi/180;
% Sim.IC = [start_slope, start_slope, 0, 0, zeros(1, GA.Sim.Con.stDim)];
wSim.Con = wSim.Con.HandleEvent(1, wSim.IC(wSim.ConCo));
if strcmp(wSim.Con.name, 'Matsuoka')
    wSim.Con = wSim.Con.Adaptation();
    
    %%%%%%%%%%%%%%% Separate Matsuoka simulation %%%%%%%%%%%%%%%%
    % Run this code to see how the neuronal controller converges
    wMOSim = deepcopy(wSim);
    
    % Initial conditions
    IC0 = wMOSim.IC(wMOSim.ConCo);
    wMOSim.Con.startup_t = 0;
%     wMOSim.Con.s_in = 3;
%     wMOSim.Con = wMOSim.Con.Adaptation(0);
    
    % tspan
    t_start = 0;
    t_end = 10;
    t_step = 0.025;
    t_span = t_start:t_step:t_end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     % TODO: don't do the MAtsuoka Sim if we only have GA!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%
	% Simulation settings:
	options = odeset('MaxStep',t_step/10,'RelTol',.5e-7,'AbsTol',.5e-8);
	
    % run simulation the first time:
    Matsuoka_tic = tic;
    [TTemp,XTemp] = ode45(@MatsDerivative,t_span,IC0,options);
    Matsuoka_runTime = toc(Matsuoka_tic);
    
    % % % use regression NN to change the gene:
    if ~isempty(GA.NN_reg) && ~isempty(GA.NN_reg_Fcn) ...
            && strcmp(wSim.Con.name, 'Matsuoka')
        % Use NN to select best value for tau gene
        thisSeq = GA.NN_reg_Fcn(GA.Gen, GA.NN_reg, thisSeq, XTemp, TTemp);
        wSim = GA.Gen.Decode(wSim,thisSeq);
        wSim.Con = wSim.Con.HandleEvent(1, wSim.IC(wSim.ConCo));
        wSim.Con = wSim.Con.Adaptation();
    end    

    % % % use classifier NN to get good genes to replace gene:
    if ~isempty(GA.NN_classi) && ~isempty(GA.NN_classi_Fcn) ...
            && strcmp(wSim.Con.name, 'Matsuoka')
        % thisSeq = GA.NN_classi_Fcn(GA.Gen, GA.NN_classi, thisSeq, XTemp, TTemp);
        
        % % % Testing New code 15/10/2017:
        lastGen = max(1,GA.Progress-1);
        lastGenes = [GA.Gen.Mutate(GA.Seqs(1:GA.Fittest(2),:,lastGen));...
            GA.Gen.Mutate(GA.Seqs(1:GA.Fittest(2),:,lastGen));...
            GA.Gen.Mutate(GA.Seqs(1:GA.Fittest(2),:,lastGen))];% add mutated copies of the last genes
        
        thisSeq = GA.NN_classi_Fcn(GA.Gen, GA.NN_classi,...
            thisSeq, lastGen,...
            lastGenes, XTemp, TTemp);
        % % % % % % % % % % % % % % % % % % %
        
        wSim = GA.Gen.Decode(wSim,thisSeq);
        wSim.Con = wSim.Con.HandleEvent(1, wSim.IC(wSim.ConCo));
        wSim.Con = wSim.Con.Adaptation();
    end    
    
    % % Time rescaling
    if ~isempty(GA.rescaleFcn)
        % Option to consider: use feedforward NN to predict the period.
        
        % Run simulation again to do the rescaling
        %   The NN might have change the period.
        [TTemp,XTemp] = ode45(@MatsDerivative,t_span,IC0,options);
        
        % use Rescaling to change the gene:
        thisSeq = GA.rescaleFcn(GA.Gen, thisSeq, XTemp, TTemp);
        wSim = GA.Gen.Decode(wSim,thisSeq);
        wSim.Con = wSim.Con.HandleEvent(1, wSim.IC(wSim.ConCo));
        wSim.Con = wSim.Con.Adaptation();
    end
    
end

    function xdot = MatsDerivative(t,x)
        xdot = wSim.Con.Derivative(t,x);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% Run the simulation
Sim_tic = tic;
wSim = wSim.Run();
Sim_runTime = toc(Sim_tic);

% check the end_Type of the simulation (why did it stop)
sim_endCond = wSim.Out.Type;

% check the ration of T_end/Sim.out.Tend:
T = wSim.Out.T;
Tend_ratio = T(end)/wSim.Out.Tend;

wSim.Con.startup_t = 0;
wSim.ConvProgr = [0, 0];
wSim.minMaxDiff = [1, 0];

% % Some more simulation initialization
% wSim.Con = wSim.Con.Reset(wSim.IC(wSim.ConCo));
% wSim.Con.FBType = 2;
% wSim.Con = wSim.Con.HandleExtFB(wSim.IC(wSim.ModCo),...
%                 wSim.IC(wSim.ConCo),wSim.Env.SurfSlope(wSim.Mod.xS));
% wSim.Env = wSim.Env.Set('Type','inc','start_slope',start_slope*180/pi);

FitInd = GA.FitFcn(:,1); % Index for fit functions
FitFcn = GA.FitFcn(:,2); % Handles for fit functions
NFit = size(GA.FitFcn,1);

% Calculate the genome's fitness
thisFit = zeros(1,max(cell2mat(FitInd')));
thisOuts = cell(1,NFit);
for f = 1:NFit
    % Preprocessing for ZMPFit
    if ~isempty(strfind(func2str(FitFcn{f}),'ZMPFit'))
        % Prepare all the required vectors
        % (torques, state, etc) and put them in wSim.Out

        % ZMP Fit should be the last one as it uses the
        % output from all previous fitness functions
        Outs = find(~cellfun(@isempty,thisOuts));
        for o = 1:length(Outs)
            wSim.Out = wSim.JoinOuts(thisOuts{Outs(o)});
        end
        wSim.Con.FBType = 2;
    end

    % Call the fitness function
    [thisFit(FitInd{f}),thisOuts{f}] = FitFcn{f}(wSim);

    % Postprocessing for VelFit
    if ~isempty(strfind(func2str(FitFcn{f}),'VelFit'))
        % Switch direction if the model walks backwards
        if thisFit(FitInd{f})<0
            revSeq = GA.Gen.SwitchDir(thisSeq);
            [Res,revSeq] = GA.Gen.CheckGenome(revSeq);
%             disp(['Switched direction ',int2str(ID)])
%             disp(GA.Gen.seq2str(thisSeq));

            if Res{1} && ~all(revSeq == 0)
                thisSeq = revSeq;
                thisFit(FitInd{f}) = -thisFit(FitInd{f});
            end
        end                    
    end

    %%%% REMOVED FOR THE TAGA-LIKE CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % TODO: return this when going back to the general case!!!
%     % Postprocessing for VelRangeFit
%     if ~isempty(strfind(func2str(FitFcn{f}),'VelRangeFit')) && ...
%         ~any(strcmp(GA.Gen.Keys(1,:),'ks_tau'))
%             % ^ This genome uses a single gain for positive and
%             % negative s_in, thus a change is irrelevant
%         % If high-level signal only works in one direction,
%         % copy those parameters for the other direction
%         newSeq = thisSeq;
%         if abs(thisFit(3))>0 && thisFit(5)==0
%             newSeq(15:18) = 0.33*newSeq(11:14);
%         end
%         if abs(thisFit(5))>0 && thisFit(3)==0
%             newSeq(11:14) = 0.33*newSeq(15:18);
%         end
%         thisSeq = newSeq;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Prepare output
switch nargout
    case 0
        if exist('ID','var')
            disp(['Genome ',num2str(ID),' results: ',num2str(thisFit)]);
        else
            disp(['Genome results: ',num2str(thisFit)]);
        end
        disp(wSim.Out.Text)
        varargout = {};
    case 2
        varargout = {thisFit, thisSeq};
	case 4
		varargout = {thisFit, thisSeq,Matsuoka_runTime,Sim_runTime};
    case 6
		varargout = {thisFit, thisSeq,Matsuoka_runTime,...
            Sim_runTime,sim_endCond,Tend_ratio};
end

end

