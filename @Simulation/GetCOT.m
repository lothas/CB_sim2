function COT = GetCOT(sim, work_in_out, count_pe)
%GETCOT Calculates and returns the Cost Of Transport for the simulation

if nargin<3
    if nargin<2
        % How to deal with positive and negative work
        work_in_out = [1 -1];
    end
    % Take into account a change in potential energy?
    count_pe = 0;
end

COT = [];
if isempty(sim.Out)
    return
end

X = sim.Out.X;
T = sim.Out.T;
Torques = sim.Out.Torques;

% Get initial and final hip position
Hip0 = zeros(2,1); Hip1 = zeros(2,1);
sim.Mod.xS = sim.Out.SuppPos(1,1);
sim.Mod.yS = sim.Out.SuppPos(1,2);
[Hip0(1), Hip0(2)] = sim.Mod.GetPos(X(1,sim.ModCo),'Hip');
sim.Mod.xS = sim.Out.SuppPos(end,1);
sim.Mod.yS = sim.Out.SuppPos(end,2);
[Hip1(1), Hip1(2)] = sim.Mod.GetPos(X(end,sim.ModCo),'Hip');

Weight = sim.Mod.GetWeight();

% Calculate distance travelled
DistanceTravelled = abs(Hip1(1)-Hip0(1));

if DistanceTravelled<3*sim.Mod.L && sim.Out.nSteps<5
    % Robot didn't walk a minimum of 3 body lengths
    return;
else
    % Calculate control effort
    StTrq = Torques(:,1)-Torques(:,2);      % Ankle torque
%     StTrq = Torques(:,1);
    StAngVel = X(:,sim.ModCo(3));           % Ankle ang. velocity
    
    SwTrq = Torques(:,2);                   % Hip torque
%     SwAngVel = X(:,sim.ModCo(4))-X(:,sim.ModCo(3));
    SwAngVel = X(:,sim.ModCo(4));           % Hip ang. velocity
    
    AnkWork = StTrq.*StAngVel;
    HipWork = SwTrq.*SwAngVel;
    
    AnkPos = AnkWork>0;
    HipPos = HipWork>0;
    
    % Calculate absolute energy input from controller
    ControlEffort = trapz(T, AnkWork.*...
        (work_in_out(1)*AnkPos + work_in_out(2)*~AnkPos)) + ...
                    trapz(T, HipWork.*...
        (work_in_out(1)*HipPos + work_in_out(2)*~HipPos));

    % Calculate Cost Of Transport
    if count_pe
        % Calculate difference in potential energy
        dPotentialE = Weight*(Hip1(2)-Hip0(2));

        COT = (ControlEffort-dPotentialE)/(Weight*DistanceTravelled);
    else
        COT = ControlEffort/(Weight*DistanceTravelled);
    end
end

