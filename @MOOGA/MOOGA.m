classdef MOOGA
    % MOOGA - Multi-Objective Optimization Genetic Algorithm
    % This class evolves controllers for the compass biped
    % model using multi-objective optimization
    
    properties
        % Optimization properties
        Generations;
        Population;
        Fittest;
        
        % Objects
        Gen;        % Genome
        Sim;        % Simulation
        Graphics = 1;   % Show/Hide simulation animation
        
        % Genomes
        Seqs;       % Genome sequences
        Fit;        % Genome fitnesses
        
        % Additional genes optimized by previous MOOGA
        BaseGen;    
        BaseSeq;
        
        % Optimization progress
        Progress = 1;   % Last generation calculated
        ReDo = 0;       % Start optimization from generation 1
        
        % Fitness functions
        NFit;
        FitFcn;     % Fitness functions handles
        
        % External function (called at the end of every generation)
        GenerationFcn;
        
        % Filenames
        FileIn;     % Input file with previous genomes (same MOOGA)
        FileOut;    % Output file
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
                                  GA.Population-2*TopPop]; % Fittest children
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
                    GA = GA.InitSim();
                    GA = GA.InitGen();
                    GA = GA.Run();
            end            
        end
        
        %% %%%%%%%%%%%%% Built-in fitness functions %%%%%%%%%%%%% %%
        function fit = HeightFit(GA,Sim) %#ok<INUSL>
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
        end
        
        function fit = VelFit(GA,Sim) %#ok<INUSL>
            X = Sim.Out.X;
            T = Sim.Out.T;
            Nt = length(T);
            L = Sim.Mod.L;
            
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
            MaxVel = L; % 1 leg length per second
            LocalMaxDist = MaxVel*T(end); % shorter than global if robot fell
            GlobalMaxDist = MaxVel*Sim.Out.Tend;

            % Give points up to local max dist
            % Normalize by global max dist (based on tend)
            if abs(MaxDist(2))>MaxDist(1)
                fit = max(MaxDist(2),-LocalMaxDist)/GlobalMaxDist;
            else
                fit = min(MaxDist(1),LocalMaxDist)/GlobalMaxDist;
            end
        end
        
        function [fit] = NrgEffFit(GA,Sim) %#ok<INUSL>
            X = Sim.Out.X;
            T = Sim.Out.T;
            Torques = Sim.Out.Torques;
            Hip0 = zeros(2,1); Hip1 = zeros(2,1);
            [Hip0(1), Hip0(2)] = Sim.Mod.GetPos(X(1,Sim.ModCo),'Hip');
            [Hip1(1), Hip1(2)] = Sim.Mod.GetPos(X(end,Sim.ModCo),'Hip');
            Weight = Sim.Mod.GetWeight();
            
            % Calculate distance travelled
            DistanceTravelled = abs(Hip1(1)-Hip0(1));
            
            if DistanceTravelled<3*Sim.Mod.L
                fit = 0;
            else
                % Calculate absolute control effort
                StTrq = Torques(:,1)-Torques(:,2);
                StAngVel = X(:,Sim.ModCo(3));
                SwTrq = Torques(:,2);
                SwAngVel = X(:,Sim.ModCo(4));
                ControlEffort = trapz(T,abs(StTrq.*StAngVel)) + ...
                                trapz(T,abs(SwTrq.*SwAngVel));

                % Calculate difference in potential energy
                dPotentialE = Weight*9.81*(Hip1(2)-Hip0(2));

                % Calculate Cost Of Transport
                COT=(ControlEffort-dPotentialE)/(Weight*DistanceTravelled);

                % Low COT (as low as 0) is best so points are given by
                fit=1/(1+5*COT);
                % COT of 0 gives 1
                % COT of 0.03 gives 0.869
                % COT of 0.12 gives 0.625
                % COT of 0.3 gives 0.4
            end
        end
    end
        
end

