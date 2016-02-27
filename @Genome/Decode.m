function [ Sim ] = Decode( Ge, Sim, Seq )
%DECODE Decode loads a string of numbers (genome) into a simulation
%   Decode receives a string of numbers and converts it based on the
%   genome keys in order to update the properties of a simulaiton
%   (including the Model, Environment, Controller, initial conditions, etc)

Res = Ge.CheckGenome(Seq);
if Res{1} == 0
%     disp(Res{2});
    % Instead of throwing an error, let's get a random sequence and
    % continue running
%     error('ERROR: Invalid sequence');
%     return;
    Seq = Ge.RandSeq();
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
                            + sign(KE(v))*Seq(SeqPos-1+abs(KE(v)));
                    end
                end
            else
                Sim.IC = Seq(SeqPos:SeqPos+Ge.Keys{2,k});
            end
            
        %% %%%%%%%%%% Controller keys %%%%%%%%%% %%
        case Sim.Con.SetKeys
            Sim.Con = Sim.Con.Set(Ge.Keys{1,k},Seq(SeqPos:SeqPos+Ge.Segments(k)-1));
        case 'Pulses'
            % For pulses Ge.Keys{2,k} is a 2 item vector
            % 1st item is the number of pulses for a joint
            % 2nd item is the joint number
            for p = 1:Ge.Keys{2,k}(1)
                P0 = SeqPos+(p-1)*Ge.KeyLength.Pulses;
                if Sim.Con.FBType == 2
                    Sim.Con = Sim.Con.AddPulse(...
                        'joint',Ge.Keys{2,k}(2),...
                        'amp',Seq(P0),...
                        'offset',Seq(P0+1),...
                        'dur',Seq(P0+2),...
                        'k_u',Seq(P0+3),...
                        'k_d',Seq(P0+4));
                else
                    Sim.Con = Sim.Con.AddPulse(...
                        'joint',Ge.Keys{2,k}(2),...
                        'amp',Seq(P0),...
                        'offset',Seq(P0+1),...
                        'dur',Seq(P0+2));
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
                        'amp',Seq(P0),...
                        'offset','ext',...
                        'dur',Seq(P0+1),...
                        'k_u',Seq(P0+2),...
                        'k_d',Seq(P0+3));
                else
                    Sim.Con = Sim.Con.AddPulse(...
                        'joint',Ge.Keys{2,k}(2),...
                        'amp',Seq(P0),...
                        'offset','ext',...
                        'dur',Seq(P0+1));
                end
            end
        case 'IC_matsuoka'
            % When this key shows up in a genome, it tells the decoder to
            % set the simulation's initial condition to a value that aids
            % the oscillators' convergence
            Sim = Sim.Init();
            wfe = mean(mean(abs(Sim.Con.win(Sim.Con.win~=0))));
%             Sim.Con.u0 = Sim.Con.beta/abs(wfe);
%             u_ss = Sim.Con.u0/Sim.Con.beta*abs(wfe);
%             u_ss = Sim.Con.Amp(1)/Sim.Con.beta*abs(wfe);
            u_ss = max(Sim.Con.Amp)/Sim.Con.beta*abs(wfe);
            MIC = zeros(1,Sim.Con.stDim);
%             MIC(1) = 10;
%             MIC(3) = u_ss;
            Sim.IC(Sim.ConCo) = MIC;
            
        
        %% %%%%%%%%%% Model keys %%%%%%%%%% %%
        case Sim.Mod.SetKeys
            Sim.Mod = Sim.Mod.Set(Ge.Keys{1,k},Seq(SeqPos:SeqPos+Ge.Segments(k)-1));
        
        %% %%%%%%%%%% Environment keys %%%%%%%%%% %%
        case Sim.Env.SetKeys
            Sim.Env = Sim.Env.Set(Ge.Keys{1,k},Seq(SeqPos:SeqPos+Ge.Segments(k)-1));
    end
    
    % Move the sequence reading position
    SeqPos = Ge.AdvSeq(SeqPos,k);
end

end

