function [ Sim ] = Decode( Ge, Sim, Genes )
%DECODE Decode loads a string of numbers (genome) into a simulation
%   Decode receives a string of numbers and converts it based on the
%   genome keys in order to update the properties of a simulaiton
%   (including the Model, Environment, Controller, initial conditions, etc)

if length(Genes)~=Ge.Length
    error(['Decode failed. Genes provided: ',num2str(length(Genes)),...
        '. Genes required: ',num2str(Ge.Length)]);
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
        
        %% %%%%%%%%%% Model keys %%%%%%%%%% %%
        
        %% %%%%%%%%%% Environment keys %%%%%%%%%% %%
    end
    
    % Move the sequence reading position
    if isfield(Ge.KeyLength,Ge.Keys{1,k})
        SeqPos = SeqPos + ...
            Ge.KeyLength.(Ge.Keys{1,k})*Ge.Keys{2,k};
    else
        SeqPos = SeqPos + Ge.Keys{2,k};
    end
end


end

