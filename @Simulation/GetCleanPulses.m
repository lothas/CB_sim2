function [varargout] = GetCleanPulses(sim)
%GETCLEANPULSES Separates the "impulse" pulses from regular pulses        
    T = sim.Out.T;
    NT = length(T);
    X = sim.Out.X;
    Torques = sim.Out.Torques;
    Impulses = zeros(size(sim.Out.Torques,1),1);
    
    if ~isempty(sim.Con.ExtPulses)
        PushOffInd = sim.Con.ExtPulses;

        % Remove "push-off" torques
        stepTime = find(diff(sim.Out.SuppPos(:,1))~=0);
        stepTime = [0; stepTime; NT];
        pulseEnd = zeros(length(stepTime),1);

        % Save lastPhi (otherwise the slope adaptation gets discarded)
        Temp = sim.Con.lastPhi;
        sim.Con.lastPhi = sim.Out.Slopes(1);

        % Maximum torque allowed as a pulse
        if ~isempty(sim.Con.MinSat)
            maxPls = max(max(abs(sim.Con.MinSat(2:end)),...
                             abs(sim.Con.MaxSat(2:end))));
        else
            Pulses = setdiff(length(sim.Con.Amp),PushOffInd);
            maxPls = abs(sim.Con.Amp(Pulses));
        end

        for s = 1:length(stepTime)-1
            thisT = T(stepTime(s)+1:stepTime(s+1));

            % Update push-off parameters
            sim.Con = sim.Con.HandleExtFB(X(stepTime(s)+1,sim.ModCo),...
                                          X(stepTime(s)+1,sim.ConCo),...
                                          sim.Out.Slopes(stepTime(s)+1));
            POdur = 0.999*sim.Con.Duration(PushOffInd)/sim.Con.omega;
            POamp = sim.Con.Amp(PushOffInd);

            if length(thisT)>1
                for ind = 1:length(PushOffInd)
                    next = find(thisT>=thisT(1)+POdur(ind),1,'first');
                    if ~isempty(next)
                        pulseEnd(s) = stepTime(s) + next;
                        POids = stepTime(s) + ...
                            find(abs(Torques(stepTime(s)+1:pulseEnd(s),1))>...
                                 abs(POamp(ind))-maxPls);
                        Torques(POids,1) = Torques(POids,1)-POamp(ind);
                        Impulses(POids) = Impulses(POids)+POamp(ind);
                    end
                end
            else
                pulseEnd(s) = stepTime(s) + 1;
                for ind = 1:length(PushOffInd)
                    POids = stepTime(s) + ...
                        find(abs(Torques(stepTime(s)+1:pulseEnd(s),1))>...
                             abs(POamp(ind))-maxPls);
                    Torques(POids,1) = Torques(POids,1)-POamp(ind);
                    Impulses(POids) = Impulses(POids)+POamp(ind);
                end
            end

        end
        sim.Con.lastPhi = Temp;
    end

    switch nargout
        case 1
            varargout{1} = Torques;
        case 2
            varargout{1} = Torques;
            varargout{2} = Impulses;
    end
end