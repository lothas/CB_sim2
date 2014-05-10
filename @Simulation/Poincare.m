function [EigVal,EigVec] = Poincare( sim )
% Calculates the Linearized Poincare eigenvalues
% Version 0.1 - 10/05/2014

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
end

% Run the simulations
PMsim = sim;
PMsim.EndCond = [1,sim.Period];
for d = Coords
    PMsim.IC = dIC(:,d);
    PMsim.Init();
    PMsim = PMsim.Run();

    dICp(:,d) = PMsim.ICstore(:,1);
end

% Calculate deviation
DP = 1/sim.PMeps*(dICp(Coords,:) - IC(Coords,:));
[EigVec,EigVal] = eig(DP,'nobalance');
EigVal = diag(EigVal);
end