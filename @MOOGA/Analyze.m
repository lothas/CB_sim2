function [ Data ] = Analyze( GA, varargin )
%ANALYZE Analyze the performance of the selected evolved controller
%   Run simulations with the selected controller over a range of slopes
%   until it is unable to walk. Calculate performance parameters such as
%   initial conditions, limit cycles, required torque/motor power,
%   required ZMP, etc.

MaxTries = 10;
base_d = 0.05;

Type = 'CL';
switch nargin 
    case 2
        Generation = GA.Progress;
        ID = varargin{1};
    case 3
        if any(isletter(varargin{2}))
            Generation = GA.Progress;
            ID = varargin{1};
            Type = varargin{2};
        else
            Generation = varargin{1};
            ID = varargin{2};
        end
    case 4
        Generation = varargin{1};
        ID = varargin{2};
        Type = varargin{3};
    otherwise
        Generation = GA.Progress;
        TopIDs = GA.GetTopPop(GA.Fittest(1));
        ID = randsample(TopIDs,1);
end

Filename = ['Gen',int2str(ID),Type,'.mat'];
if exist(Filename,'file') == 2
    In = load(Filename);
    Data = In.Data;
    if Data.Done
        switch Type
            case {'OL','CL'}
                DispRes();
            case 's_in'
        end
        return
    end
    
    switch Type
        case 's_in'
            start_slope = 0;
            start_s_in = Data.s_in(end);
            if length(Data.s_in)<2
                ds_in = base_d;
            else
                ds_in = sign(Data.s_in(end))*median(abs(diff(Data.s_in)));
            end
        case {'OL','CL'}
            start_slope = Data.Slopes(end);
            if length(Data.Slopes)<2
                d_slope = base_d;
            else
                d_slope = sign(Data.Slopes(end))*median(abs(diff(Data.Slopes)));
            end
    end
else
    % Initialize output
    switch Type
        case 's_in'
            Data.s_in = []; % High level speed command
        case {'OL','CL'}
            Data.Slopes = []; % Slopes
    end
    Data.IC = []; % Initial conditions
    Data.LCx = {}; % Limit cycle states
    Data.LCt = {}; % Limit cycle time
    Data.LCtorques = {}; % Limit cycle Torques
    Data.MTorques = []; % Max torques
    Data.Power = []; % Actuator power required
    Data.EigV = []; % Eigenvalues
    Data.Period = []; % Period
    Data.StepLength = [];
    Data.Speed = [];
    Data.LCGRF = []; % Limit cycle GRF
    Data.LCZMP = []; % Limit cycle ZMP
    Data.MuFric = []; % Max friction coefficient
    Data.MZMP = []; % Max ZMP
    
    % Energies:
    Data.dPotE = []; % Potential energy gained/lost
    Data.dKinEIm = []; % Kinetic energy lost to impact
    Data.dKinEFr = []; % Kinetic energy lost to friction
    Data.ContEff = []; % Control effort
    Data.AbsContEff = []; % Absolute control effort
    Data.COT = []; % Cost of transport
    
    start_slope = 0;
    d_slope = 0;
    start_s_in = 0;
    ds_in = 0;
    
    Data.Done = 0;
    
    % Store genome information
    Data.Gen = GA.Gen;
    Data.Seq = GA.Seqs(ID,:,Generation);    
end

% Run simulation starting from slope 0 and then increasing/decreasing slope
Sim = deepcopy(GA.Sim);
Sim.EndZMP = 1; % Stop for ZMP (1) - don't stop (0)
Sim.Graphics = 0;
if Sim.Graphics == 1
    Sim.Fig = figure();
end
Sim.PMFull = 1; % Run poincare map on all 5 coords
switch Type
    case {'OL','s_in'}
        Sim.Con.FBType = 0;
    case 'CL'
        Sim.Con.FBType = 2;
end

Sim = Data.Gen.Decode(Sim, Data.Seq);

disp('Performing analysis...')
while ~Data.Done
    Tries = 0;
    while Tries<MaxTries
        % Run simulation
        switch Type
            case {'OL','CL'}
                Asim = WalkOnSlope(Sim,Data,start_slope,...
                                            start_slope+d_slope);
            case 's_in'
                Asim = WalkAtSpeed(Sim,Data,start_s_in,...
                                            start_s_in+ds_in);
        end

        if Asim.Out.Type == 5
            % Simulation converged
            fprintf(' - ');
            cprintf('*green','OK!\n')
            
            switch Type
                case {'OL','CL'}
                    Data = SavePerformance(Asim,start_slope+d_slope,Data);

                    % Keep going to next slope
                    if d_slope == 0
                        d_slope = base_d;
                    else
                        start_slope = start_slope+d_slope;
                        d_slope = d_slope/abs(d_slope)*base_d;
                    end
                case 's_in'
                    Data = SavePerformance(Asim,start_s_in+ds_in,Data);

                    % Keep going to next speed
                    if ds_in == 0
                        ds_in = base_d;
                    else
                        start_s_in = start_s_in+ds_in;
                        ds_in = ds_in/abs(ds_in)*base_d;
                    end
            end
            break
        else
            fprintf(' - ');
            cprintf('*red','FAILED!\n')
            disp(Asim.Out.Text);
            
            switch Type
                case {'OL','CL'}
                    if Asim.Out.Type ~= 0
                        % Try running from 0 to the required slope
                        Asim = WalkOnSlope(Sim,Data,0,start_slope+d_slope);

                        if Asim.Out.Type == 5
                            % Simulation converged
                            fprintf(' - ');
                            cprintf('*green','OK!\n')

                            Data = SavePerformance(Asim,start_slope+d_slope,Data);

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
                case 's_in'
                    % Keep going to next speed
                    if ds_in>0
                        ds_in = ds_in+base_d;
                    else
                        ds_in = ds_in-base_d;
                    end
                    Tries = Tries + 1;
            end     
        end
    end
    
    if Tries == MaxTries
        switch Type
            case {'OL','CL'}
                % Simulation failed to converge
                % Start checking negative slopes
                if d_slope > 0
                    start_slope = 0;
                    d_slope = -base_d;
                else
                    Data.Done = 1;
                    break;
                end
            case 's_in'
                % Simulation failed to converge
                % Start checking negative speeds
                if ds_in > 0
                    start_s_in = 0;
                    ds_in = -base_d;
                else
                    Data.Done = 1;
                    break;
                end
        end
    end
end

% Sort data
switch Type
    case {'OL','CL'}
        Last = find(diff(Data.Slopes)<0,1,'first');
        NS = length(Data.Slopes);
    case 's_in'
        Last = find(diff(Data.s_in)<0,1,'first');
        NS = length(Data.s_in);
end

Fields = fieldnames(Data);
for f = 1:length(Fields)
    if strcmp(Fields{f},'Done')
        continue
    end
    
    [r,c] = size(Data.(Fields{f}));
    if c == NS
        Data.(Fields{f}) = [fliplr(Data.(Fields{f})(:,Last+1:end)),...
                            Data.(Fields{f})(:,1:Last)];
    else
        if r == NS
            Data.(Fields{f}) = [flipud(Data.(Fields{f})(Last+1:end,:));...
                                Data.(Fields{f})(1:Last,:)];
        else
            disp('Error: wrong number of data points');
        end
    end
end

switch Type
    case {'OL','CL'}
        % Calculate Zones
        Zones = {};
        % Period 1 zones
        dslope = mean(diff(Data.Slopes));
        splits = find(diff(Data.Slopes)>1.0001*dslope);
        if isempty(splits)
            Zones{1} = {1:length(Data.Slopes)};
        else
            P1zone = {1:splits(1)};
            for zn = 1:length(splits)-1
                P1zone(1,zn+1) = {splits(zn)+1:splits(zn+1)};
            end
            P1zone(1,end+1) = {splits(end)+1:length(Data.Slopes)};
            Zones{1} = P1zone;
        end
        % Period 2 zones
        Zones{2} = {};
        for z1 = 1:length(Zones{1})
            P2IDs = intersect(find(Data.Period(:,1)==2),Zones{1}{z1});
            if isempty(P2IDs)
                continue
            end
            splits = find(diff(P2IDs)>1);
            if isempty(splits)
                Zones{2} = [Zones{2}, P2IDs(1:length(P2IDs))];
            else
                P1zone = {P2IDs(1:splits(1))};
                for zn = 1:length(splits)-1
                    P1zone(1,zn+1) = {P2IDs(splits(zn)+1:splits(zn+1))};
                end
                P1zone(1,end+1) = {P2IDs(splits(end)+1:length(P2IDs))}; %#ok<AGROW>
                Zones{2} = [Zones{2}, P1zone];
            end
        end
        if isempty(Zones{2})
            Zones(2) = [];
        end
        Data.Zones = Zones;
    case 's_in'
        Data.Zones = 1;
end


% Sort eigenvalues lines
% [Data.EigV] = SortLines(Data.EigV);
[Data.EigV(1:5,:)] = ConnectVectors(Data.EigV(1:5,:));
if length(Data.Zones)>1
    [Data.EigV(6:10,:)] = ConnectVectors(Data.EigV(6:10,:));
end
                
save(Filename,'Data');
switch Type
    case {'OL','CL'}
        DispRes();
    case 's_in'
end

    function sim = WalkOnSlope(sim_in, Data, start_s, end_s)
        sim = deepcopy(sim_in);
        sim  = sim.Init();
        fprintf('Processing data on %.2f degrees slope',end_s);
        
        if start_s == 0
            if strcmp(sim.Con.name, 'Matsuoka')
                sim.IClimCyc = 0*sim.IC;
                slope = sim.Env.start_slope;
                sim.IClimCyc(sim.ModCo(1:2)) = [slope, -slope];
                sim.IClimCyc(sim.ConCo(1)) = 1;
            else
                sim.IClimCyc = 0*sim.IC;
            end
        else
            sim.IClimCyc = Data.IC(1:sim.stDim,Data.Slopes == start_s);
        end
        
        sim = sim.WalkOnSlope(start_s, end_s, ~start_s*5, 75);
    end

    function sim = WalkAtSpeed(sim_in, Data, start_s_in, s_in)
        sim = deepcopy(sim_in);
        sim  = sim.Init();
        fprintf('Processing data for s_in = %.2f', s_in);
        
        if s_in == 0
            sim.IClimCyc = 0*sim.IC;
        else
            sim.IClimCyc = Data.IC(1:sim.stDim,Data.s_in == start_s_in);
        end
        
        sim = sim.WalkAtSpeed(s_in, 75);
    end

    function Data = SavePerformance(sim,S,Data)
        switch Type
            case {'OL','CL'}
                Data.Slopes(end+1) = S;
                Slope = S;
            case 's_in'
                Data.s_in(end+1) = S;
                Slope = 0;
                sim.Con.s_in = S;
        end
        Data.Period(end+1,1) = sim.Period(1);
        % Initial conditions (for each step until a full period)
        Coords = 1:sim.stDim*sim.Period(1);
        Data.IC(Coords,end+1) = ...
            reshape(sim.ICstore(:,1:sim.Period(1)),[],1); 
        
        % Calculate eigenvalues
        if sim.PMFull
            EachP = sim.stDim;
        else
            EachP = length(sim.ModCo);
        end
        EigV = zeros(EachP*sim.Period(1),1);
        for p = 1:sim.Period(1)
            sim.IClimCyc = sim.ICstore(:,p);
            Coords = 1+EachP*(p-1):EachP*p;
            [EigV(Coords),~] = sim.Poincare();
        end
        Data.EigV(1:length(EigV),end+1) = EigV;
        
        % Prepare simulation for single period evaluation
        sim = sim.SetTime(0,0.003,5);
        sim.Mod.xS = 0; sim.Mod.yS = 0;
        sim.Env = sim.Env.Set('Type','inc','start_slope',Slope);
        sim.EndCond = [1,sim.Period(1)]; % Run for one full period
        sim = sim.Init();
        sim.Mod.LegShift = sim.Mod.Clearance;
        sim.IC = Data.IC(1:sim.stDim,end); % sim.IClimCyc;
        sim.Con = sim.Con.HandleExtFB(sim.IC(sim.ModCo),...
                sim.IC(sim.ConCo),sim.Env.SurfSlope(sim.Mod.xS));
        sim.Con = sim.Con.Reset(sim.IC(sim.ConCo));
    
        sim = sim.Run();
        
        Data.LCx{end+1} = sim.Out.X; % Limit cycle states
        Data.LCt{end+1} = sim.Out.T; % Limit cycle time
        Data.LCtorques{end+1} = sim.Out.Torques; % Limit cycle torques
        
        Data.Period(end,2) = max(sim.Out.T); % Step duration
        X0 = Data.IC(1:sim.stDim,end);
        Data.StepLength(end+1) = norm(sim.Mod.GetPos(X0,'NS') - ...
                                      sim.Mod.GetPos(X0,'S'));
        Data.Speed(end+1) = Data.StepLength(end)/Data.Period(end,2);
        
        % Deal with torques
        % Separate impulses
        [Torques, Impulses] = sim.GetCleanPulses();
        Torques = [Torques, Impulses];
        % Calculate max torque (positive or negative)
        mT = min(Torques,[],1); MT = max(Torques,[],1);
        Moverm = find(abs(MT)>=abs(mT));
        moverM = find(abs(MT)<abs(mT));
        Data.MTorques(end+1,Moverm) = MT(Moverm);
        Data.MTorques(end,moverM) = mT(moverM);
        
        % Actuator power required
        AnkPower = Torques(:,1).*sim.Out.X(:,3);
        HipPower = Torques(:,2).*(sim.Out.X(:,4)-sim.Out.X(:,3));
        Data.Power(end+1,:) = [max(abs(AnkPower)); max(abs(HipPower))];
        
        % ZMP
        NT = length(sim.Out.T);
        % Process the data
        GRF = zeros(NT,2);
        ZMP = zeros(NT,1);
        for t=1:NT
            thisX = sim.Out.X(t,sim.ModCo);
            sim.Mod.Torques = sim.Out.Torques(t,:)';
            GRF(t,:) = sim.Mod.GetGRF(thisX)';
            ZMP(t) = Torques(t,1)/GRF(t,2); % Ankle torque/GRFy
        end

        Data.LCGRF{end+1} = GRF;
        Data.LCZMP{end+1} = ZMP;
        % Max coefficient of friction required and minimum GRFy
        Data.MuFric(end+1,:) = [max(abs(GRF(:,1)./GRF(:,2))); min(GRF(:,2))];
        Data.MZMP(end+1,:) = [max(ZMP); min(ZMP)]; % ZMP front; back
        
        % Calculate energies
        Data = CalcEnergy(sim,Data);
        
        % Save progress
        save(Filename,'Data');
    end

    function DispRes()
        disp(['Min/Max slope: ',...
        num2str(min(Data.Slopes),'%.2f'),' / '...
        num2str(max(Data.Slopes),'%.2f')])
    end

    function [Lines] = ConnectVectors( Lines )
        N=size(Lines,2);

        for d=1:N-1
            L=length(find(Lines(:,d)~=0));
            if L==10
                % To shorten the computation time of calculating
                % all the permutations of a 10 elements vector,
                % we remove the 2 eigenvalues that are closest
                % to zero (there's a zero eigenvalue in the PM)

                CurVec=[Lines(1:L,d), (1:L)']; % Vector that holds the current points
                NextVec=[Lines(1:L,d+1), (1:L)']; % Vector that holds the next set of points

                CurVecS=sortrows(CurVec,1);
                NextVecS=sortrows(NextVec,1);

                % Grab the last 8
                CurVec=CurVecS(3:end,1);
                NextVec=NextVecS(3:end,1);

                L=8;
                if length(CurVec)~=length(NextVec)
                    A=1;
                end
            else
                CurVec=Lines(1:L,d); % Vector that holds the current points
                NextVec=Lines(1:L,d+1); % Vector that holds the next set of points
            end

            % Calculate all possible permutations of the next points
            % This takes a LONG LONG LONG while when M is large
            PermVec=perms(NextVec);

            % Calculate distance from current vector to all permutations
            Dist=abs(PermVec-repmat(CurVec',size(PermVec,1),1));

            % Find the most distant point for each permutation
            MaxDist=max(Dist,[],2);

            % Find the permutation where the max. distance is the smallest
            MinPerm=find(MaxDist==min(MaxDist));

            MP=length(MinPerm);
            if MP==1
                NextVec=PermVec(MinPerm,:)';
                if L~=8
                    Lines(1:L,d+1)=NextVec;
                else
                    Lines(CurVecS(1:2,2),d+1)=Lines(NextVecS(1:2,2),d+1);
                    Lines(CurVecS(3:end,2),d+1)=NextVec;
                end
            else            
                % Too many best permutations found
                % We'll select the one with the minimal overall distance
                BestPerms=PermVec(MinPerm,:);

                % Calculate absoulte distance from CurVec to each candidate
                Dist=abs(BestPerms-repmat(CurVec',MP,1));
                AbsDist=diag(Dist*Dist');

                MinPerm=find(AbsDist==min(AbsDist));
                if length(MinPerm)>1
                    % screw it, choose the first one
                    NextVec=BestPerms(MinPerm(1),:)';
                else
                    NextVec=BestPerms(MinPerm,:)';
                end

                if L~=8
                    Lines(1:L,d+1)=NextVec;
                else
                    Lines(CurVecS(1:2,2),d+1)=Lines(NextVecS(1:2,2),d+1);
                    Lines(CurVecS(3:end,2),d+1)=NextVec;
                end
            end
        end
    end

    function [Lines] = SortLines( Lines )
       % Works great as long as lines don't intersect
       N=size(Lines,2);

       % We'll only look at lines 1 and 6 and switch the whole sets
       % 1:5 and 6:10 accordingly. Let row one be always "above" row 6
       for cl=1:N
           if Lines(1,cl)<Lines(6,cl)
               % Need to switch
               Temp=Lines(1:5,cl:end);
               Lines(1:5,cl:end)=Lines(6:10,cl:end);
               Lines(6:10,cl:end)=Temp;
           end
       end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Potential, Kinetic, and other energies calculation %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function Data = CalcEnergy(Sim,Data)
        ESim = deepcopy(Sim);
        
        % Get limit cycle
        X = ESim.Out.X;
        T = ESim.Out.T;
        Torques = ESim.Out.Torques;
        
        % Get initial and final position
        Hip0 = zeros(2,1); Hip1 = zeros(2,1);
        
        ESim.Mod.xS = ESim.Out.SuppPos(1,1);
        ESim.Mod.yS = ESim.Out.SuppPos(1,2);
        [Hip0(1), Hip0(2)] = ESim.Mod.GetPos(X(1,ESim.ModCo),'Hip');
        
        ESim.Mod.xS = ESim.Out.SuppPos(end,1);
        ESim.Mod.yS = ESim.Out.SuppPos(end,2);
        [Hip1(1), Hip1(2)] = ESim.Mod.GetPos(X(end,ESim.ModCo),'Hip');
        
        Weight = ESim.Mod.GetWeight(); % Get weight

        % Calculate distance travelled
        DistanceTravelled = abs(Hip1(1)-Hip0(1));

        % Calculate absolute control effort
        StTrq = Torques(:,1)-Torques(:,2);
        StAngVel = X(:,ESim.ModCo(3));
        SwTrq = Torques(:,2);
        SwAngVel = X(:,ESim.ModCo(4));
        AbsControlEffort = trapz(T,abs(StTrq.*StAngVel)) + ...
                           trapz(T,abs(SwTrq.*SwAngVel));
        ControlEffort = trapz(T,StTrq.*StAngVel) + ...
                        trapz(T,SwTrq.*SwAngVel);

        % Calculate difference in potential energy
        dPotentialE = Weight*(Hip1(2)-Hip0(2));

        % Calculate kinetic energy lost to impact
        if ESim.Period(1)>1
            warning('Period doubling calculation not implemented')
            dKineticEIm = 0;
        else
            Xb = X(end,ESim.ModCo);
            Xa = ESim.Mod.CalcImpact(Xb);
            if norm(Xa-X(1,ESim.ModCo))>1e-6
                warning('Erroneous limit cycle')
            end
            K0 = ESim.Mod.GetKineticEnergy(Xb);
            K1 = ESim.Mod.GetKineticEnergy(Xa);
            dKineticEIm = K1 - K0;
        end
        
        % Calculate kinetic energy lost to friction
        StFric = - ESim.Mod.dampA*StAngVel...
                 - ESim.Mod.dampH*(StAngVel-SwAngVel);
        SwFric = - ESim.Mod.dampH*(SwAngVel-StAngVel);
        dKineticEFr = trapz(T,StFric.*StAngVel) + ...
                      trapz(T,SwFric.*SwAngVel);
        
        % Calculate Cost Of Transport
        COT=(ControlEffort-dPotentialE)/(Weight*DistanceTravelled);
        
        % Save data
        Data.dPotE(end+1,:) = dPotentialE;
        Data.dKinEIm(end+1,:) = dKineticEIm;
        Data.dKinEFr(end+1,:) = dKineticEFr;
        Data.ContEff(end+1,:) = ControlEffort;
        Data.AbsContEff(end+1,:) = AbsControlEffort;
        Data.COT(end+1,:) = COT;
    end

end

