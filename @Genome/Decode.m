function [ Sim ] = Decode( Ge, Sim, Genes )
%DECODE Decode loads a string of numbers (genome) into a simulation
%   Decode receives a string of numbers and converts it based on the
%   genome keys in order to update the properties of a simulaiton
%   (including the Model, Environment, Controller, initial conditions, etc)

if ~Ge.CheckGenome(Genes)
    error('ERROR: Invalid sequence');
%     return;
end

SeqPos = 1; % Position along the genome sequence
for k = 1:size(Ge.Keys,2)
    % Go over each key
    switch Ge.Keys{1,k}
        %% %%%%%%%%%% Simulation keys in general %%%%%%%%%% %%
        case {'IC','ic','init cond','initial conditions'}
            if isfield(Ge.KeyExtra,Ge.Keys{1,k})
                KE = Ge.KeyExtra.(Ge.Keys{1,k});
                % Store initial conditions provided based
                % on the KeyExtra value for 'IC'
                for v = 1:length(KE)
                    if KE(v) ~= 0
                        Sim.IC(v) = Sim.IC(v) ...
                            + sign(KE(v))*Genes(SeqPos-1+abs(KE(v)));
                    end
                end
            else
                Sim.IC = Genes(SeqPos:SeqPos+Ge.Keys{2,k});
            end
            
        %% %%%%%%%%%% Controller keys %%%%%%%%%% %%
        case Sim.Con.SetKeys
            Sim.Con = Sim.Con.Set(Ge.Keys{1,k},Genes(SeqPos:SeqPos+Ge.Keys{2,k}(1)-1));
        case 'Pulses'
            % For pulses Ge.Keys{2,k} is a 2 item vector
            % 1st item is the number of pulses for a joint
            % 2nd item is the joint number
            for p = 1:Ge.Keys{2,k}(1)
                P0 = SeqPos+(p-1)*Ge.KeyLength.Pulses;
                if Sim.Con.FBType == 2
                    Sim.Con = Sim.Con.AddPulse(...
                        'joint',Ge.Keys{2,k}(2),...
                        'amp',Genes(P0),...
                        'offset',Genes(P0+1),...
                        'dur',Genes(P0+2),...
                        'k_u',Genes(P0+3),...
                        'k_d',Genes(P0+4));
                else
                    Sim.Con = Sim.Con.AddPulse(...
                        'joint',Ge.Keys{2,k}(2),...
                        'amp',Genes(P0),...
                        'offset',Genes(P0+1),...
                        'dur',Genes(P0+2));
                end
            end
        case 'ExtPulses'
            % For pulses Ge.Keys{2,k} is a 2 item vector
            % 1st item is the joint number
            % 2nd item is the number of pulses for that joint
            for p = 1:Ge.Keys{2,k}(1)
                P0 = SeqPos+(p-1)*Ge.KeyLength.ExtPulses;
                if Sim.Con.FBType == 2
                    Sim.Con = Sim.Con.AddPulse(...
                        'joint',Ge.Keys{2,k}(2),...
                        'amp',Genes(P0),...
                        'offset','ext',...
                        'dur',Genes(P0+1),...
                        'k_u',Genes(P0+2),...
                        'k_d',Genes(P0+3));
                else
                    Sim.Con = Sim.Con.AddPulse(...
                        'joint',Ge.Keys{2,k}(2),...
                        'amp',Genes(P0),...
                        'offset','ext',...
                        'dur',Genes(P0+1));
                end
            end
        
        %% %%%%%%%%%% Model keys %%%%%%%%%% %%
        case Sim.Mod.SetKeys
            Sim.Mod = Sim.Mod.Set(Ge.Keys{1,k},Genes(SeqPos:SeqPos+Ge.Keys{2,k}(1)-1));
        
        %% %%%%%%%%%% Environment keys %%%%%%%%%% %%
        case Sim.Env.SetKeys
            Sim.Env = Sim.Env.Set(Ge.Keys{1,k},Genes(SeqPos:SeqPos+Ge.Keys{2,k}(1)-1));
    end
    
    % Move the sequence reading position
    SeqPos = Ge.AdvSeq(SeqPos,k);
end


end

