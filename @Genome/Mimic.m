function [ Ge ] = Mimic( Ge, Seq1, Seq2, Sim, Slope )
%MIMIC Calculates the gains needed to mimic Ge2 behavior on Slope
%   Mimic updates the controller's gains to mimic the torque signals
%   of Ge2 at a certain slope.
%   NOTE: Slope should be provided in degrees.

Sim1 = deepcopy(Sim);
Sim2 = deepcopy(Sim);

Sim1 = Ge.Decode(Sim1, Seq1);
Sim2 = Ge.Decode(Sim2, Seq2);

K = Sim1.Con.MimicGains(Sim2.Con,Slope);

% Re-enconde the genome with the new gains
end

