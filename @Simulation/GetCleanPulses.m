function [varargout] = GetCleanPulses(sim, varargin)
%GETCLEANPULSES Separates the "impulse" pulses from regular pulses   
    if nargin == 1
        T = sim.Out.T;
        X = sim.Out.X;
        Torques = sim.Out.Torques;
    else
        T = varargin{1};
        X = varargin{2};
        Torques = varargin{3};
    end
    
    Impulses = zeros(size(Torques,1),1);
    if ~isempty(sim.Con.ExtPulses)
        ImpInd = sim.Con.ExtPulses;
        [ImpJoint,~,~] = find(sim.Con.OutM(:,sim.Con.ExtPulses));
        Impulses = zeros(size(Torques,1),length(ImpJoint));

        % Remove "push-off" torques
        NT = length(T);
        stepTime = find(diff(sim.Out.SuppPos(:,1))~=0);
        stepTime = [0; stepTime; NT];
        pulseEnd = zeros(length(stepTime),1);

        % Save lastPhi (otherwise the slope adaptation gets discarded)
        Temp = sim.Con.lastPhi;
        sim.Con.lastPhi = sim.Out.Slopes(1);

        % Maximum torque allowed as a pulse
        maxPls = zeros(length(ImpInd),1);
        if ~isempty(sim.Con.MinSat)
            for p = 1:length(ImpInd)
                % Find pulses applied to the same joint
                Pulses = find(sim.Con.OutM(ImpJoint(p),:)==1);
                Pulses = setdiff(Pulses,ImpInd(p));
                maxPls(p) = max(max(abs(sim.Con.MinSat(Pulses)),...
                                    abs(sim.Con.MinSat(Pulses))));
            end
        else
            for p = 1:length(ImpInd)
                % Find pulses applied to the same joint
                Pulses = find(sim.Con.OutM(ImpJoint(p),:)==1);
                Pulses = setdiff(Pulses,ImpInd(p));
                if ~isempty(Pulses)
                    maxPls(p) = abs(sim.Con.Amp(Pulses));
                else
                    maxPls(p) = 0;
                end
            end
        end

        for s = 1:length(stepTime)-1
            thisT = T(stepTime(s)+1:stepTime(s+1));

            % Update push-off parameters
            sim.Con = sim.Con.HandleExtFB(X(stepTime(s)+1,sim.ModCo),...
                                          X(stepTime(s)+1,sim.ConCo),...
                                          sim.Out.Slopes(stepTime(s)+1));
            POdur = 0.999*sim.Con.Duration(ImpInd)/sim.Con.omega;
            POamp = sim.Con.Amp(ImpInd);

            if length(thisT)>1
                for ind = 1:length(ImpInd)
                    next = find(thisT>=thisT(1)+POdur(ind),1,'first');
                    if ~isempty(next)
                        pulseEnd(s) = stepTime(s) + next;
                        POids = stepTime(s) + ...
                            find(abs(Torques(stepTime(s)+1:pulseEnd(s),...
                                             ImpJoint(ind)))>=...
                                     0.999*abs(POamp(ind))-maxPls);
                        Torques(POids,ImpJoint(ind)) = ...
                            Torques(POids,ImpJoint(ind))-POamp(ind);
                        Impulses(POids,ImpJoint(ind)) = ...
                            Impulses(POids,ImpJoint(ind))+POamp(ind);
                    end
                end
            else
                pulseEnd(s) = stepTime(s) + 1;
                for ind = 1:length(ImpInd)
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