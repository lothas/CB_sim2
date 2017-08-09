classdef MOOGA
    % MOOGA - Multi-Objective Optimization Genetic Algorithm
    % This class evolves controllers for the compass biped
    % model using multi-objective optimization
    
    properties
        % Optimization properties
        Generations;
        Population;
        Fittest;
        
        WeightedPairing = 1; % Individuals have higher chances of
                             % "mating" when they belong to higher
                             % pareto fronts
        JOAT = 1;   % "Jack of all trades"
                    % Set to 0 uses the normal pareto front approach
                    % Set to 1 takes the first pareto front from the
                    %   general population and the remaining fronts from
                    %   a quantile of the population that is fit in
                    %   many aspects.
                    % Set to 2 takes all fronts from the quantile of the
                    %   population, a.k.a. "master of none".
        Quant = 0.8;    % Drops 20% of the population that are bad at one
                        % or more fitness aspects (for JOAT>0)

        % Use NN to modify some genes, akin to fancy epigenetic gene
        % expression. How cool is that?
        NN = [];
        NNFcn = [];
        
        % Rescale some genes to get a certain temporal performance
        rescaleFcn = [];
        
        % Objects
        Gen;        % Genome
        Sim;        % Simulation
        Graphics = 1;   % Show/Hide simulation animation
        
        % Genomes
        Seqs;       % Genome sequences
        Fit;        % Genome fitnesses
        Parents;    % Parent IDs for each genome
        
        % Additional genes optimized by previous MOOGA
        BaseGen;
        BaseSeq;
        
        % Optimization progress
        Progress = 0;   % Last generation calculated
        CompTime = 0;   % Computation time elapsed
        ReDo = 0;       % Start optimization from generation 1
        
        % Fitness functions
        NFit;
        FitFcn;     % Fitness functions handles
        FitIDs;     % Which values to use from within Fit
        FitMinMax;  % Which fit values to show/find as larger/smaller than
                
        % External function (called at the end of every generation)
        GenerationFcn;
        
        % Filenames:
        FileIn;     % Input file with previous genomes (same MOOGA)
        FileOut;    % Output file
		
		% Run times:
        totGenTime = [] % save the generation time
		MLseqRunTime = [];	% run time for the Matsuoka CPG by itself
		simRunTime = [];	% run time for the CPG with the CB
        
        % some more parameters for the paper:
        sim_endCond = [];
        Tend_ratio = [];
		
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function GA = MOOGA(varargin)
            GA.Gen.COType = '2point';
            GA.Gen.MutType = 'norm';
            switch nargin
                case 2
                    GA.Generations = varargin{1};
                    GA.Population = varargin{2};
                    TopPop = floor(GA.Population/5);
                    GA.Fittest = [TopPop,... % Fittest copy
                                  TopPop,... % Mutated fittest copy
                                  ceil(0.025*GA.Population)]; % New random genomes
                case 3
                    GA.Generations = varargin{1};
                    GA.Population = varargin{2};
                    GA.Fittest = varargin{3};
                otherwise
                    GA.Generations = 5;
                    GA.Population = 1000;
                    TopPop = floor(GA.Population/5);
                    GA.Fittest = [TopPop,... % Fittest copy
                                  TopPop,... % Mutated fittest copy
                                  GA.Population-2*TopPop]; % Fittest children
                    Keys = {'Height','Weight','Age';
                                   1,       1,    1};
                    Range = {0, 0, 0;
                             1, 1, 1};
                    GA.Gen = Genome(Keys, Range);
                    GA.GenerationFcn = @GA.Test;
                    GA = GA.InitGen();
                    GA = GA.Run();
            end            
        end
        
        function GA = SetFittest(GA,Top,MutTop,NewRnd)
            GA.Fittest = [floor(Top/100*GA.Population),...
                          floor(MutTop/100*GA.Population),...
                          floor(NewRnd/100*GA.Population)];
        end
        
        function varargaout = Find(GA,varargin)
            Gnrtn = min(GA.Progress+1,GA.Generations);
            switch nargin
                case 1
                    find_max = true; %#ok<NASGU>
                case 2
                    Reqs = varargin{1};
                case 3
                    Reqs = varargin{1};
                    Gnrtn = varargin{2};
            end
            
            PFits = repmat(GA.FitMinMax, ...
                           size(GA.Fit(:,:,Gnrtn), 1), 1) ...
                    .*GA.Fit(:,:,Gnrtn);
            MaxFits = GA.FitMinMax.*max(PFits);
            
            if exist('find_max','var')
                Max = cell(2,max(cell2mat(GA.FitFcn(:,1)')));
                for f=1:size(GA.FitFcn,1)
                    FitInd = GA.FitFcn{f,1};
                    Max{1,FitInd(1)} = ...
                        MOOGA.GetFitFcnName(GA.FitFcn{f,2});
                    Max(2,FitInd) = ...
                        num2cell(MaxFits(FitInd));
                end
                disp(Max)
                return
            end
                    
            [R,~] = size(Reqs);
            if R>1
                % Only some conditions provided as pairs:
                % [Fitness number, Minimum required]
                Reqs = zeros(max(cell2mat(GA.FitFcn(:,1)')),1);
                for r = 1:R
                    Reqs(varargin{1}(r,1)) = varargin{1}(r,2);
                end
            else
                if length(Reqs)~=size(GA.FitFcn,1)
                    disp('Number of fitness values is incorrect');
                    return;
                end
            end

            % Find results that fit the requirements
            Conds = ones(GA.Population,1);
            for f = 1:max(cell2mat(GA.FitFcn(:,1)'))
                if sign(GA.FitMinMax(f)) >= 0
                    Conds = Conds & ...
                        GA.Fit(:,f,Gnrtn)>=Reqs(f);
                else
                    Conds = Conds & ...
                        GA.Fit(:,f,Gnrtn)<=-Reqs(f);
                end 
            end

            Fits = find(Conds);
            if length(Fits)>15 && nargout<1
                disp([int2str(length(Fits)),...
                    ' results fit the requirements']);
                return
            else
                if isempty(Fits)
                    disp('No matches found');
                    return
                end
                out = [Fits, GA.Fit(Fits,:,Gnrtn)];
                switch nargout
                    case 0
                        disp(num2str(out,'%.3f'));
                    case 1
                        varargaout = out;
                end
            end
        end
        
        function ShowFinalDist(GA,only_best,gen)
            if nargin<2
                only_best = 1;
            end
            if nargin<3
                gen = GA.Progress;
            end
            
            if only_best == 1
                top_ids = GA.GetTopPop(GA.Fittest(1));
                data = GA.Seqs(top_ids,:,gen);
            else
                data = GA.Seqs(:,:,gen);
            end
            
            Nfeat = size(data,2);
            cols = ceil(sqrt(Nfeat));
            rows = ceil(Nfeat/cols);
            Nsamp = size(data,2); Nbins = max(10,ceil(Nsamp/20));
            for i = 1:Nfeat
                subplot(rows,cols,i);
                hist(data(:,i),Nbins);
            end                
        end
        
        function varargout = PrintRuntime(GA)         
            if GA.CompTime > 0
                Time = GA.CompTime;
                hours = floor(Time/3600); Time = Time - hours*3600;
                minutes = floor(Time/60); Time = Time - minutes*60;
                seconds = Time;
                % Prepare output
                outStr = 'Total runtime: ';
                if hours > 1
                    outStr = [outStr, int2str(hours),' hours, '];
                else
                    if hours == 1
                        outStr = [outStr, '1 hour, '];
                    end
                end
                if minutes > 1
                    outStr = [outStr, int2str(minutes),' minutes, '];
                else
                    if minutes == 1
                        outStr = [outStr, '1 minute, '];
                    end
                end
                if seconds == 1
                    outStr = [outStr, '1 second'];               
                else
                    if seconds > 0
                        outStr = [outStr, num2str(seconds),' seconds'];
                    else
                        % no seconds, remove ', '
                        outStr = outStr(1:end-2);
                    end
                end
            else
                outStr = 'Total runtime: unknown';
            end
            
            if nargout == 0
                disp(outStr)
            else
                varargout = {outStr};
            end
        end
    end
    
    methods(Static)
        %% %%%%%%%%%%%%% Built-in fitness functions %%%%%%%%%%%%% %%
        function [fit,out] = HeightFit(Sim)
            X = Sim.Out.X;
            T = Sim.Out.T;
            Nt = length(T);
            L = Sim.Mod.L;
            Points = zeros(1,Nt);
            for t = 1:Nt
                % Get hip height at each instant
                [~, Y] = Sim.Mod.GetPos(X(t,Sim.ModCo),'Hip');
                
                % Award or substract points based on the hip height
                % For the simple model the hip can't be higher than leg length
                % So we'll focus on height above/below a certain threshold

                % Max points between 100% and 85% of leg height L
                % Points go from 1 to 0 between 85% and 70%
                % (85% gives an aperture between the legs of 63 deg)
                % (70% gives an aperture between the legs of 90 deg)
                Points(t) = max(min(1,(Y-0.7*L)/(0.15*L)),0);
            end
            fit = trapz(T,Points);
            % Normalize fit
            fit = min(fit/trapz(linspace(0,Sim.Out.Tend-2*Sim.tstep_normal,Nt),...
                            ones(size(T))),1);
            out = [];
        end
        
        function [fit,out] = VelFit(Sim)
            X = Sim.Out.X;
            T = Sim.Out.T;
            Nt = length(T);
            L = Sim.Mod.L;
            Sim.Mod.xS = Sim.Out.SuppPos(1,1);
            Sim.Mod.yS = Sim.Out.SuppPos(1,2);
            
            % Find the farthest point that the robot reached
            MaxDist = [0,0];
            for t = 1:Nt
                [StanceX,~] = Sim.Mod.GetPos(X(t,Sim.ModCo),'S');
                [SwingX,~] = Sim.Mod.GetPos(X(t,Sim.ModCo),'NS');
                Fwd = Sim.Out.SuppPos(t,1)+min(StanceX,SwingX);
                Bwd = Sim.Out.SuppPos(t,1)+max(StanceX,SwingX);
                MaxDist = [max(MaxDist(1),Fwd), min(MaxDist(2),Bwd)];
            end
            
            % We want to encourage the robot to move forward
            % But we don't want it to move too fast
            % So we'll limit the points it can get based on a max velocity
            MaxVel = 1.5*L; % 1.5 leg length per second
            LocalMaxDist = MaxVel*T(end); % shorter than global if robot fell
            GlobalMaxDist = MaxVel*Sim.Out.Tend;

            % Give points up to local max dist
            % Normalize by global max dist (based on tend)
            if abs(MaxDist(2))>MaxDist(1)
                fit = max(MaxDist(2),-LocalMaxDist)/GlobalMaxDist;
            else
                fit = min(MaxDist(1),LocalMaxDist)/GlobalMaxDist;
            end
            out = Sim.Out;
        end
        
        function [fit,out] = NrgEffFit(Sim)
            % Calculate Cost Of Transport
            COT = Sim.GetCOT([1,0], 0);
            
            if isempty(COT)
                fit = 0;
            else
                % Low COT (as low as 0) is best so points are given by
                fit=1/(1+5*COT);
            end
            % COT of 0 gives 1
            % COT of 0.03 gives 0.869
            % COT of 0.12 gives 0.625
            % COT of 0.3 gives 0.4
            out = [];
        end
        
        function [fit,out] = EigenFit(Sim)
            if isempty(Sim.Period)
                fit = 0;
            else
                % Calculate numerical poincare
                [EigVal,~] = Sim.Poincare();
                fit = 1-max(abs(EigVal));
                if fit<0
                    fit = 0;
                end
            end
            out = [];
        end
        
        function [Sim] = SetSimSlope(Sim,vfit,UD)
            parK = min(max(0.005/vfit,0.003),0.01);
            leadway = 1;

            Sim.IC = Sim.ICstore(:,1);
            Sim.Env = ...
                Sim.Env.Set('Type','inf','start_slope',0,...
                                 'parK',UD*parK,'start_x',leadway);
            Sim.Con.FBType = 2;
            Sim.Mod.LegShift = Sim.Mod.Clearance;
            Sim = Sim.SetTime(0,0.1,180);
            Sim.Con = Sim.Con.Reset(Sim.IC(Sim.ConCo));
            Sim.Con = Sim.Con.HandleExtFB(Sim.IC(Sim.ModCo),...
                Sim.IC(Sim.ConCo),Sim.Env.SurfSlope(Sim.Mod.xS));
            Sim = Sim.Init();
        end
        
        function [fit,out] = UphillFitRun(Sim)
            SlopeSim = deepcopy(Sim);
            
            DistanceTravelled = Sim.Out.SuppPos(end,1);
            SlopeSim.Mod.xS = 0; SlopeSim.Mod.yS = 0;
            
            if DistanceTravelled<3*SlopeSim.Mod.L
                fit = 0;
            else
                V = DistanceTravelled/SlopeSim.Out.T(end);
                SlopeSim = MOOGA.SetSimSlope(SlopeSim,V,1);
                SlopeSim = SlopeSim.Run();
                fit = MOOGA.UphillFit(SlopeSim);
            end
            out = SlopeSim.Out;
        end
        
        function [fit,out] = UphillFit(Sim)
            fit = Sim.MaxSlope;
            out = [];
        end
        
        function [fit,out] = DownhillFitRun(Sim)
            SlopeSim = deepcopy(Sim);
            
            DistanceTravelled = Sim.Out.SuppPos(end,1);
            SlopeSim.Mod.xS = 0; SlopeSim.Mod.yS = 0;
            
            if DistanceTravelled<3*SlopeSim.Mod.L
                fit = 0;
            else
                V = DistanceTravelled/SlopeSim.Out.T(end);
                SlopeSim = MOOGA.SetSimSlope(SlopeSim,V,-1);
                SlopeSim.Con = SlopeSim.Con.Reset(SlopeSim.IC(SlopeSim.ConCo));
                SlopeSim = SlopeSim.Run();
                fit = MOOGA.DownhillFit(SlopeSim);
            end
            out = SlopeSim.Out;
        end
                
        function [fit,out] = DownhillFit(Sim)
            fit = -Sim.MinSlope;
            out = [];
        end
        
        function [fit,out] = ZMPFit(Sim)
            if Sim.Out.SuppPos(end,1)<3*Sim.Mod.L
                fit = 0;
                out = [];
            else
                X = Sim.Out.X;
                Torques = Sim.GetCleanPulses();
                
                [ZMPfront,ZMPback] = Sim.Mod.GetZMP(X, Torques);
                
                % Max foot size
                MaxFront = 0.25; % cm (ankle to toes)
                MaxBack = 0.15; % cm (ankle to heel)
                if abs(ZMPfront)/MaxFront>abs(ZMPback)/MaxBack
                    Max = MaxFront;
                    ZMP = abs(ZMPfront);
                else
                    Max = MaxBack;
                    ZMP = abs(ZMPback);
                end
                % Make the fitness function nonlinear:
                % Grade changes slowly around 0, fast around Max
                % Max grade is 1 at 0, 0.5 at Max, goes to 0 in infinity.
                fit = 1./(1+(ZMP/Max).^3);
                out = [];
            end
        end
        
        function [fit,out] = SlopeFit(Sim,dir)
            % Start running up/downwards at intervals and checks that the
            % simulation converges on each slope before moving forward
            
            Slope = 0;
            dSlope = dir*1;
            SlSim = deepcopy(Sim);
            SlSim.doGoNoGo = 2;
            SlSim.GNGThresh = [5,10];
            if isempty(SlSim.IClimCyc)
                SlSim.IClimCyc = 0*SlSim.IC;
                Leadway = 3;
            else
                Leadway = 0;
            end
            out = [];
            while 1
                SlSim = SlSim.WalkOnSlope(Slope,Slope+dSlope,Leadway,40);
                if Slope~=0
                    % Save part of the output
                    i = find(dir*SlSim.Out.Slopes*58>=dir*(Slope+dSlope),...
                        1,'first');
                else
                    % Save the whole output
                    i = length(SlSim.Out.Slopes);
                end
                out = SlSim.JoinOuts(out,i);
                
                if SlSim.Out.Type == 6
                   % Simulation GO
                   Slope = Slope+dSlope;
                else
                    % Slope = max(abs(SlSim.Out.Slopes))*58;
                    % Add the percentage of total time it managed to walk
                    % on the last slope before falling
                    Slope = Slope + SlSim.Out.T(end)/40*dSlope;
                    break;
                end
            end
            fit = min(abs(Slope),10);
        end
        
        function [fit,out] = UpSlopeFit(Sim)
            [fit,out] = MOOGA.SlopeFit(Sim,1);
        end
        
        function [fit,out] = DownSlopeFit(Sim)
            [fit,out] = MOOGA.SlopeFit(Sim,-1);
        end
        
        function [fit,out] = UDZMPFit(Sim,dir)
            % Calls the up/down fitness function and then calculates the
            % ZMP of that run
            [Slope,out] = MOOGA.SlopeFit(Sim,dir);
            
            ZMPSim = deepcopy(Sim);
            ZMPSim.Out = out;
            [ZMP,~] = MOOGA.ZMPFit(ZMPSim);
            
%             fit = [50*ZMP*(1-cosd(Slope)), Slope, ZMP];
            fit = [(Slope+5.0*ZMP)/15.0 Slope, ZMP];
        end
        
        function [fit,out] = ZMPUpFit(Sim)
            [fit,out] = MOOGA.UDZMPFit(Sim,1);
        end
        
        function [fit,out] = ZMPDownFit(Sim)
            [fit,out] = MOOGA.UDZMPFit(Sim,-1);
        end
        
        function [fit,out] = UpDownFit(Sim)
            % Combines the upslope and downslope fitness results
            
            if isempty(Sim.Period)
                % Weed out results that didn't converge
                fit = [0 0 0];
                out = [];
                return;
            end
                
            % Run the upslope test
            [up_fit,up_out] = MOOGA.SlopeFit(Sim,1);
            % Run the downslope test
            [down_fit,down_out] = MOOGA.SlopeFit(Sim,-1);
            
            fit = [up_fit*down_fit up_fit down_fit];
            
            UDSim = deepcopy(Sim);
            UDSim.Out = up_out;
            out = UDSim.JoinOuts(down_out);
        end
        
        function [max_vel, COT, max_s, out] = DeltaVelFit(Sim, dir, base_vel)
            % Starts walking faster/slower at intervals and checks that the
            % simulation converges at each speed before moving forward
            
            max_vel = base_vel;
            max_s = 0;
            s_in = 0;
            ds_in = dir*0.5;
            VelSim = deepcopy(Sim);
            VelSim.doGoNoGo = 2;
            VelSim.GNGThresh = [5,10];
            VelSim.MinMaxStore = zeros(length(VelSim.stepsSS), ...
                                       max(VelSim.GNGThresh)+1);
            out = [];
            
            maxSingleT = 20; % Time alloted for a single run
            maxAllT = 20*maxSingleT; % Time alloted for all runs
            elapsedT = 0;    % Time elapsed in all runs
            COT = [];
            
            while 1
                % Run a simulation of the robot walking with different speed
                s_in = s_in + ds_in;
                if abs(s_in) > 20
                    % Stop increasing after 20 resoultion steps
                    break;
                end
                
                VelSim = VelSim.WalkAtSpeed(s_in, maxSingleT);
                elapsedT = elapsedT + VelSim.Out.T(end) - VelSim.Out.T(1);
                if elapsedT > maxAllT
                    break
                end
                
                out = VelSim.JoinOuts(out);
                
                if any(VelSim.Out.Type == [0,5,6])
                    % Time ran out, sim converged or was converging
                    % Simulation GO
                    % Get average velocity
                    per = Sim.GetPeriod(0.6);
                    T = per(2);
                    n = min(VelSim.Out.nSteps, size(VelSim.ICstore,2));
                    
                    if n > 5
                        % At least 5 steps should be taken
                        % (though go-no-go should take care of it as well)
                        avg_vel = sum(2*sin(abs(VelSim.ICstore(1,2:n))) ...
                            *Sim.Mod.L/T) / (n-1);
                    else
                        avg_vel = max_vel;
                    end
                    
%                     if dir*avg_vel > (1+0.005*dir)*dir*max_vel
                        % Expect at least a 0.1% increase
                    d_vel = dir*(avg_vel - max_vel);
                    if d_vel >= 0.01 && d_vel < 0.125
                        % Expect at least a 2% change in the range [0-0.5]
                        % and at most 25%
                        max_vel = avg_vel;
                        max_s = s_in;
                    else
                        break;
                    end
                else
                    break;
                end
                
                % Calculate the COT of this run
                this_COT = VelSim.GetCOT([1,0], 0);

                if isempty(COT)
                    if isempty(this_COT)
                        COT = 1; % The cost of dragging an object with coefficient
                                 % of friction = 1
                    else
                        COT = this_COT;
                    end
                else
                    if isempty(this_COT)
                        if VelSim.Out.nSteps<3
                            % Sim failed quickly on new speed, don't update
                        else
                            COT = this_COT;
                        end
                    else
                        COT = this_COT;
                    end
                end   
            end
            
            if isempty(COT)
                COT = 1; % The cost of dragging an object with coefficient
                         % of friction = 1
            end
        end
        
        function [fit,out] = VelRangeFit(Sim)
            % Combines the faster and slower velocity fitness results
            
            if isempty(Sim.Period) && Sim.StepsTaken < 10
                % Weed out results that didn't converge
                fit = [0 0 0 0 0 0 0 0];
                out = [];
                return;
            end
            
            % Get last speed achieved with s_in=0
            % Get average velocity
            per = Sim.GetPeriod(1);
            T = per(2);
            n = min(Sim.Out.nSteps, size(Sim.ICstore,2));
            base_vel = sum(2*sin(Sim.ICstore(1,2:n))*Sim.Mod.L/T) / (n-1);
            
            % Run the increased velocity test (s_in>0)
            [fast_vel, fast_COT, fast_s, fast_out] = ...
                MOOGA.DeltaVelFit(Sim, 1, base_vel);
            % Run the decreased velocity test (s_in<0)
            [slow_vel, slow_COT, slow_s, slow_out] = ...
                MOOGA.DeltaVelFit(Sim, -1, base_vel);
            
            fit = [abs(fast_vel - slow_vel), 1/(1+fast_COT+slow_COT), ...
                slow_s, fast_s, slow_vel, fast_vel, slow_COT, fast_COT];
            
            % Combine simulation outputs (in case other fitness function
            % needs them)
            FSSim = deepcopy(Sim);
            FSSim.Out = fast_out;
            out = FSSim.JoinOuts(slow_out);
        end
        
        function [name] = GetFitFcnName(fcn_handle)
            fstr = strsplit(func2str(fcn_handle),{'.','('});
            if length(fstr)==2
                name = fstr{2};
            else
                name = fstr{3};
            end
        end
    end
        
end

