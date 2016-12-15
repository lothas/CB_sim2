function [EigVal,EigVec] = Poincare( sim )
% Calculates the Linearized Poincare eigenvalues
% Version 0.1 - 10/05/2014

if isempty(sim.IClimCyc)
    EigVal = [];
    EigVec = [];
    return
end

TempG = sim.Graphics;
sim.Graphics = 0;

if sim.PMFull == 1
    Ncoord = sim.stDim;
else
    Ncoord = length(sim.ModCo);
end
Coords = 1:Ncoord;

% Limit cycle initial conditions
IC = repmat(sim.IClimCyc, 1, Ncoord);

% Disturbed initial conditions
dIC = IC;
dICp = IC;
for d = Coords
    dIC(d,d) = dIC(d,d) + sim.PMeps;
    if d == sim.ConCo
        % Check CPG after disturbance
        if dIC(d,d) > sim.Con.P_th
            dIC(d,d) = dIC(d,d) - (sim.Con.P_th - sim.Con.P_reset);
        end
    end
end


% Run the simulations
PMsim = copy(sim);
Period = sim.GetPeriod();
PMsim.EndCond = [1,Period(1)];
PMsim = PMsim.SetTime(0,2*Period(1)*Period(2));
Slope = PMsim.Env.SurfSlope(PMsim.Mod.xS);
for d = Coords
    PMsim.IC = dIC(:,d);
    PMsim = PMsim.Init();
    PMsim.Con = PMsim.Con.Reset(PMsim.IC(PMsim.ConCo));
    PMsim.Con = PMsim.Con.HandleExtFB(PMsim.IC(PMsim.ModCo)',...
        PMsim.IC(PMsim.ConCo),Slope);
    PMsim = PMsim.Run();

    if PMsim.Out.Type ~= 4
        % Robot didn't complete the step
        EigVal = 2*ones(Ncoord,1);
        EigVec = eye(Ncoord);
        return;
    end
    dICp(:,d) = PMsim.ICstore(:,1);
end

% Calculate deviation
if sim.PMFull == 1
    DeltaM = [dICp(sim.ModCo,:) - IC(sim.ModCo,:);
              sim.Con.PhaseDiff(dICp(sim.ConCo,:),IC(sim.ConCo,:))];
else
    DeltaM = dICp(Coords,:) - IC(Coords,:);
end
          
DP = 1/sim.PMeps*(DeltaM);
[EigVec,EigVal] = eig(DP,'nobalance');
EigVal = diag(EigVal);

sim.Graphics = TempG;
end