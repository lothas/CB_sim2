function [ Out ] = Analyze( GA, varargin )
%ANALYZE Analyze the performance of the selected evolved controller
%   Run simulations with the selected controller over a range of slopes
%   until it is unable to walk. Calculate performance parameters such as
%   initial conditions, limit cycles, required torque/motor power,
%   required ZMP, etc.

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
    Data = load(Filename);
    Out = Data.Out;
    start_slope = Out.Slopes(end);
    d_slope = Out.Slopes(end)-Out.Slopes(end-1);
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
    % Run simulation
    Asim = WalkOnSlope(Sim,Out,start_slope,start_slope+d_slope);
        
    if Asim.Out.Type == 5
        % Simulation converged
        Out = SavePerformance(Asim,start_slope+d_slope,Out);
        
        % Keep going to next slope
        if d_slope == 0
            d_slope = 1;
        else
            start_slope = start_slope+d_slope;
        end
    else
%         if Asim.Out.Type == 0
%             % Simulation timed-out before converging
%             % Keep going to next slope
%             if d_slope == 0
%                 d_slope = 2;
%             else
%                 start_slope = start_slope+d_slope;
%             end
%         end
        % Simulation failed to converge, start checking negative slopes
        if d_slope == 1
            start_slope = 0;
            d_slope = -1;
        else
            break;
        end
    end
end

save(Filename,'Out');

    function sim = WalkOnSlope(sim_in, Data, start_s, end_s)
        sim = deepcopy(sim_in);
        disp(['Processing data on ',num2str(end_s),' degrees slope']);
        % Simulation parameters
        sim = sim.SetTime(0,0.03,200);
        sim.EndCond = 2; % Run until converge

        % Some more simulation initialization
        sim.Mod.LegShift = sim.Mod.Clearance;
        sim.Con.FBType = 2;
        sim.Con = sim.Con.Reset();
        if start_s == 0
            sim.IC = 0*sim.IC;
            sim.Con = sim.Con.HandleEvent(1, sim.IC(sim.ConCo));
            leadway = 5;
        else
            sim.IC = Data.IC(:,Data.Slopes == start_s);
            sim.Con = sim.Con.HandleExtFB(sim.IC(sim.ModCo),sim.IC(sim.ConCo));
            leadway = 1;
        end
        sim.Env = sim.Env.Set('Type','finite','start_slope',start_s,...
                              'start_x',leadway,'end_slope',end_s);
        sim = sim.Init();

        % Simulate
        sim = sim.Run();
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
        sim.Con = sim.Con.Reset();
        sim.IC = sim.IClimCyc;
        sim.Con = sim.Con.HandleExtFB(sim.IC(sim.ModCo),sim.IC(sim.ConCo));
    
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

