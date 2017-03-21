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
if sim.PMstrob == 1
    PMsim.EndCond = 0;
    PMsim = PMsim.SetTime(0,sim.Period(1)*sim.Period(2));
else
    PMsim.EndCond = [1,sim.Period(1)];
    PMsim = PMsim.SetTime(0,2*sim.Period(1)*sim.Period(2));
end
Slope = PMsim.Env.SurfSlope(PMsim.Mod.xS);
for d = Coords
    PMsim.IC = dIC(:,d);
    PMsim = PMsim.Init();
    PMsim.Con = PMsim.Con.Reset(PMsim.IC(PMsim.ConCo));
    PMsim.Con = PMsim.Con.HandleExtFB(PMsim.IC(PMsim.ModCo)',...
        PMsim.IC(PMsim.ConCo),Slope);
    PMsim = PMsim.Run();

    if (PMsim.PMstrob == 0 && PMsim.Out.Type ~= 4) ...
            || (PMsim.PMstrob == 1 && PMsim.Out.Type ~= 0)
        % Robot didn't complete the step
        EigVal = 2*ones(Ncoord,1);
        EigVec = eye(Ncoord);
        return;
    end
    
    if PMsim.PMstrob == 0
        dICp(:,d) = PMsim.ICstore(:,1);
    else
        dICp(:,d) = PMsim.Out.X(end,:)';
    end
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