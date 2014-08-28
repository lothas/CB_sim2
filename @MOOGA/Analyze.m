function [ Out ] = Analyze( GA, varargin )
%ANALYZE Analyze the performance of the selected evolved controller
%   Run simulations with the selected controller over a range of slopes
%   until it is unable to walk. Calculate performance parameters such as
%   initial conditions, limit cycles, required torque/motor power,
%   required ZMP, etc.

MaxTries = 15;
base_d = 0.5;

switch nargin 
    case 2
        Generation = GA.Progress;
        ID = varargin{1};
    case 3
        Generation = varargin{1};
        ID = varargin{2};
    otherwise
        Generation = GA.Progress;
        TopIDs = GA.GetTopPop(GA.Fittest(1));
        ID = randsample(TopIDs,1);
end

Filename = ['Gen',int2str(ID),'.mat'];
if exist(Filename,'file') == 2
    In = load(Filename);
    Out = In.Data;
    start_slope = Out.Slopes(end);
    if length(Out.Slopes)<2
        d_slope = base_d;
    else
        d_slope = sign(Out.Slopes(end))*median(abs(diff(Out.Slopes)));
    end
else
    % Initialize output
    Out.Slopes = []; % Slopes
    Out.IC = []; % Initial conditions
    Out.LCx = {}; % Limit cycle states
    Out.LCt = {}; % Limit cycle time
    Out.LCtorques = {}; % Limit cycle Torques
    Out.MTorques = []; % Max torques
    Out.Power = []; % Actuator power required
    Out.EigV = []; % Eigenvalues
    Out.Period = []; % Period
    Out.LCZMP = []; % Limit cycle ZMP
    Out.MZMP = []; % Max ZMP
    start_slope = 0;
    d_slope = 0;
end

% Run simulation starting from slope 0 and then increasing/decreasing slope
Sim = deepcopy(GA.Sim);
Sim.Graphics = 0;
if Sim.Graphics == 1
    Sim.Fig = figure();
end
Sim.PMFull = 1; % Run poincare map on all 5 coords

Sim = GA.Gen.Decode(Sim, GA.Seqs(ID,:,Generation));

while 1
    Tries = 0;
    while Tries<MaxTries
        % Run simulation
        Asim = WalkOnSlope(Sim,Out,start_slope,start_slope+d_slope);

        if Asim.Out.Type == 5
            % Simulation converged
            fprintf(' - ');
            cprintf('*green','OK!\n')
            
            Out = SavePerformance(Asim,start_slope+d_slope,Out);

            % Keep going to next slope
            if d_slope == 0
                d_slope = base_d;
            else
                start_slope = start_slope+d_slope;
                d_slope = d_slope/abs(d_slope)*base_d;
            end
            break
        else
            fprintf(' - ');
            cprintf('*red','FAILED!\n')
            disp(Asim.Out.Text);
            
            if Asim.Out.Type ~= 0
                % Try running from 0 to the required slope
                Asim = WalkOnSlope(Sim,Out,0,start_slope+d_slope);

                if Asim.Out.Type == 5
                    % Simulation converged
                    fprintf(' - ');
                    cprintf('*green','OK!\n')

                    Out = SavePerformance(Asim,start_slope+d_slope,Out);

                    % Keep going to next slope
                    if d_slope == 0
                        d_slope = base_d;
                    else
                        start_slope = start_slope+d_slope;
                        d_slope = d_slope/abs(d_slope)*base_d;
                    end
                    break
                else
                    % Simulation failed to converge
                    fprintf(' - ');
                    cprintf('*red','FAILED!\n')
                    disp(Asim.Out.Text);

                    % Try the next slope
                    if d_slope>0
                        d_slope = d_slope+base_d;
                    else
                        d_slope = d_slope-base_d;
                    end
                    Tries = Tries + 1;
                end
            else
                % Try the next slope
                if d_slope>0
                    d_slope = d_slope+base_d;
                else
                    d_slope = d_slope-base_d;
                end
                Tries = Tries + 1;
            end
        end
    end
    
    if Tries == MaxTries
        % Simulation failed to converge, start checking negative slopes
        if d_slope > 0
            start_slope = 0;
            d_slope = -base_d;
        else
            break;
        end
    end
end

% Sort data
Data = Out;
LastUp = find(diff(Data.Slopes)<0,1,'first');
NSlopes = length(Data.Slopes);
Fields = fieldnames(Data);
for f = 1:length(Fields)
    [r,c] = size(Data.(Fields{f}));
    if c == NSlopes
        Data.(Fields{f}) = [fliplr(Data.(Fields{f})(:,LastUp+1:end)),...
                            Data.(Fields{f})(:,1:LastUp)];
    else
        if r == NSlopes
            Data.(Fields{f}) = [flipud(Data.(Fields{f})(LastUp+1:end,:));...
                                Data.(Fields{f})(1:LastUp,:)];
        else
            disp('Error: wrong number of data points');
        end
    end
end

    function sim = WalkOnSlope(sim_in, Data, start_s, end_s)
        sim = deepcopy(sim_in);
        fprintf('Processing data on %.2f degrees slope',end_s);
        
        if start_s == 0
            sim.IClimCyc = 0*sim.IC;
        else
            sim.IClimCyc = Data.IC(:,Data.Slopes == start_s);
        end
        
        sim = sim.WalkOnSlope(start_s, end_s, ~start_s*5, 75);
    end

    function Data = SavePerformance(sim,Slope,Data)
        Data.Slopes(end+1) = Slope;
        Data.IC(:,end+1) = sim.IClimCyc; % Initial conditions
        Data.Period(end+1,1) = sim.Period(1);
        
        % Calculate eigenvalues
        [Data.EigV(:,end+1),~] = sim.Poincare();
        
        % Prepare simulation for single step evaluation
        sim = sim.SetTime(0,0.003,5);
        sim.Env = sim.Env.Set('Type','inc','start_slope',Slope);
        sim.EndCond = [1,sim.Period(1)]; % Run until converge
        sim = sim.Init();
        sim.Mod.LegShift = sim.Mod.Clearance;
        sim.IC = sim.IClimCyc;
        sim.Con = sim.Con.Reset(sim.IC(sim.ConCo));
        sim.Con = sim.Con.HandleExtFB(sim.IC(sim.ModCo),...
                sim.IC(sim.ConCo),sim.Env.SurfSlope(sim.Mod.xS));
    
        sim = sim.Run();
        
        Data.LCx{end+1} = sim.Out.X; % Limit cycle states
        Data.LCt{end+1} = sim.Out.T; % Limit cycle time
        Data.Period(end,2) = max(sim.Out.T);
        Data.LCtorques{end+1} = sim.Out.Torques; % Limit cycle torques
        % Separate impulses
        [Torques, Impulses] = sim.GetCleanPulses();
        ImpTorqs = [Impulses, Torques];
        % Calculate max torque (positive or negative)
        mT = min(ImpTorqs,[],1); MT = max(ImpTorqs,[],1);
        Moverm = find(abs(MT)>=abs(mT));
        moverM = find(abs(MT)<abs(mT));
        Data.MTorques(end+1,Moverm) = MT(Moverm);
        Data.MTorques(end,moverM) = mT(moverM);
        
        % Actuator power required
        AnkPower = ImpTorqs(:,2).*sim.Out.X(:,3);
        HipPower = ImpTorqs(:,3).*(sim.Out.X(:,4)-sim.Out.X(:,3));
        Data.Power(end+1,:) = [max(abs(AnkPower)); max(abs(HipPower))];
        
        % ZMP
        NT = length(sim.Out.T);
        % Process the data
        GRF = zeros(NT,2);
        ZMP = zeros(NT,1);
        for t=1:NT
            thisX = sim.Out.X(t,sim.ModCo);
            GRF(t,:) = sim.Mod.GetGRF(thisX)';
            ZMP(t) = Torques(t,1)/GRF(t,2); % Ankle torque/GRFy
        end

        Data.LCZMP{end+1} = ZMP;
        Data.MZMP(end+1,:) = [max(ZMP); min(ZMP)]; % ZMP front; back
        
        % Save progress
        save(Filename,'Data');
    end
end

