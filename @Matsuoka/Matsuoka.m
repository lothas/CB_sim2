classdef Matsuoka < handle & matlab.mixin.Copyable
    % Version 0.1 - 04/01/2016
    
    % Matsuoka oscillator class
    
    properties
        name = 'Matsuoka'
        
        % Parameters
        tau = 1;
        tav = 1;
        beta = 0.1;
        u0 = 1;
        win = [0, 3; 1, 0];
        wex = [];
        W = [];
        
        stDim = 4; % state dimension
        nEvents = 1; % num. of simulation events
        
        % Controller Output
        startup_t = 0;
        nPulses = 1;
        OutM = [0, 0; 1 -0.1];
        Amp0 = 10;
        Amp = 10;
        ExtPulses = [];
        
        % Saturation
        MinSat; MaxSat;
        
        % Feedback
        FBType = 0;  % 0 - no feedback
        
        s_in = 0;    % Speed input - higher level input to control the
                     % desired walking speed
        
        % Gains
        ks_tau = 0;
        ks_out = 0;
        
        % Phase reset
        ExtP_reset = []; % set to a certain phase to use phase reset
                         % on foot contact
        
        % Set keys
        SetKeys = {'npulses', 'nneurons', 'n_pulses', 'n_neurons', ...
            'tau', 'tau_u', 'tav', 'tau_v', 'beta', 'win', 'wex', ...
            'amp0', 'amp', 'fbtype', 'feedback', 'fb', ...
            'ks_tau', 'speed_tau', 'tau_speed_gain', ...
            'ks_out', 'speed_out', 'torque_speed_gain'};
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function [MO] = Matsuoka(varargin)
            
%             MO.nEvents = ;
        end
        
        function MO = Reset(MO,Phase)
%             if nargin<2
%                 MO.Switch = 0*MO.Switch;
%             else
%                 % Find which Torques should be active
%                 Start = MO.Offset;
%                 End = Start + MO.Duration;
%                 [Xperc] = MO.GetPhasePerc(Phase);
%                 On = Xperc >= Start & Xperc < End;
%                 MO.Switch = (MO.Amp.*On)';
%             end 
        end
        
        function MO = SetOutMatrix(MO, ppj)
            % ppj - Pulses per joint
            % Provided as, [# of E/F neuron pairs for joint 1, ...
            %               # of E/F neuron pairs for joint 2, ... ]
            n_joints = length(ppj);
            MO.nPulses = sum(ppj);
            MO.OutM = zeros(n_joints, 2*MO.nPulses);
            p = 1;
            for j = 1:n_joints
                MO.OutM(j, p:p + 2*ppj(j) - 1) = ...
                    repmat([1, -1], 1, ppj(j));
                p = p + 2*ppj(j);
            end
        end
            
        function [Torques] = Output(MO, t, MOX, ~)
            y = max(MOX(1:2*MO.nPulses,:),0);
%             y = diag(MO.Amp)*max(MOX(1:2*MO.nPulses,:),0);
            try
                Torques = MO.OutM*y;
            catch
                disp('fuck!')
            end
%             Torques = MO.OutM*MO.Switch;

            % Apply saturation
            if ~isempty(MO.MinSat)
                for j = 1:length(MO.MinSat)
                    Torques(j,:) = ...
                        min(max(Torques(j,:),MO.MinSat(j)),MO.MaxSat(j));
                end
            end
            
            Torques(:,t<MO.startup_t) = 0*Torques(:,t<MO.startup_t);
        end

        function [per] = GetPeriod(MO)
            per = 0.754;
%             per = (MO.P_th-MO.P_reset)/MO.omega;
        end
        
        function [Xperc] = GetPhasePerc(MO,X)
            phase = atan2(X(1), X(MO.nPulses*2+1));
            Xperc = phase/2/pi;
        end
        
        function diff = PhaseDiff(MO,ph1,ph2)
            diff = ph1 - ph2;
%             % wrap it up
%             Range = MO.P_th - MO.P_reset;
%             diff(diff>Range/2) = diff(diff>Range/2) - Range;
%             diff(diff<-Range/2) = diff(diff<-Range/2) + Range;
        end
        
        function [MO] = Adaptation(MO, ~)
            if MO.FBType == 0
                % NO FEEDBACK
%                 MO.omega = MO.omega0;
%                 MO.Amp = MO.Amp0;
            else
%                 MO.omega = MO.omega0 + ...
%                     min(0,Phi)*MO.kOmega_d + ...    % Phi<0
%                     max(0,Phi)*MO.kOmega_u;         % Phi>0
%                 MO.Amp = MO.Amp0 + ...
%                     min(0,Phi)*MO.kTorques_d + ...  % Phi<0
%                     max(0,Phi)*MO.kTorques_u;       % Phi>0
            end
            
            % Apply higher-level speed input
            MO.tau = MO.tau + MO.ks_tau*MO.s_in;
            MO.tav = MO.tav + MO.ks_tau*MO.s_in;
            MO.Amp = MO.Amp0 + MO.ks_out*MO.s_in;
            MO.W = (diag(1./MO.Amp)*(diag(MO.Amp)*(MO.win+MO.wex))')';
        end
        
        % %%%%%% % Derivative % %%%%%% %
        function [Xdot] = Derivative(MO, ~, X)
            % X = [u_i; v_i];
            u = X(1:2*MO.nPulses,:);
            v = X(2*MO.nPulses+1:end);
            y = max(u,0);
            
            udot = 1/MO.tau*(-u - MO.beta*v + MO.Amp + ...
                             MO.W*y);
            vdot = 1/MO.tav*(-v+y);
            
            Xdot = [udot;
                    vdot];
        end
        
        % %%%%%% % Events % %%%%%% %
        function [value, isterminal, direction] = Events(MO, X)
            value = ones(MO.nEvents,1);
            isterminal = ones(MO.nEvents,1);
            direction = -ones(MO.nEvents,1);
            
            % Check for firing neuron
            Torques = MO.Output(0, X, 0);
            value(1) = -Torques(2,1);
% 
%             % Check for leg extension signal (for clearance)
% %             value(2) = MO.P_LegE - X;
%             
%             % Check for switching on/off signal
%             Xperc = MO.GetPhasePerc(X);
%             value(3:2+MO.nPulses) = MO.Offset - Xperc;
%             value(3+MO.nPulses:2+2*MO.nPulses) = ...
%                 MO.Offset + MO.Duration - Xperc;
        end
        
        function [MO,Xa] = HandleEvent(MO, EvID, Xb)
            Xa = Xb;
%             switch EvID
%                 case 1
%                     % Neuron fired
%                     for i=1:MO.nPulses
%                         if MO.Offset(i)==0;
%                             % Turn on signal now
%                             MO.Switch(i) = MO.Amp(i);
%                         end
%                     end
%                     Xa = MO.P_reset; % reset phase
%                 case 2
%                     % Extend the leg
%                     % This is done from the simulation file
%                 case num2cell(3:2+MO.nPulses)
%                     % Switch on signal
%                     MO.Switch(EvID-2) = MO.Amp(EvID-2);
%                 case num2cell(3+MO.nPulses:2+2*MO.nPulses)
%                     % Switch off signal
%                     PulseID = EvID-(2+MO.nPulses);
%                     MO.Switch(PulseID) = 0;
%                     if any(PulseID == MO.ExtPulses)
%                         % Set offset back to 200%
%                         MO.Offset(PulseID) = 2;
%                     end
%             end
        end
        
        function [MO, Xmod, Xcon] = HandleExtFB(MO, Xmod, Xcon, Slope)
            % This function is called when the leg hits the ground
%             if MO.FBType > 0
%                 % Perform adaptation based on terrain slope
                MO = MO.Adaptation(Slope);
%             end
%             
%             if ~isempty(MO.ExtP_reset)
%                 % Perform a phase reset
%                 Xcon = MO.ExtP_reset;
%                 
%                 % Check if any event happens at ExtP_reset
%                 [value, it, dir] = MO.Events(Xcon); %#ok<NASGU,ASGLU>
%                 EvIDs = find(value == 0);
%                 for ev = 1:length(EvIDs)
%                     [MO,Xcon] = MO.HandleEvent(EvIDs(ev),Xcon);
%                 end
%             end
%             
%             switch MO.FBImpulse
%                 case 1 % 1 - add a constant value
%                     Xmod(3:4) = Xmod(3:4) + MO.AngVelImp;
%                 case 2 % 2 - set ang. vel. to certain value
% %                     delta = MO.AngVelImp - Xmod(3:4)
%                     Xmod(3:4) = MO.AngVelImp;
%             end
%             
%             % Activate external pulses
%             MO.Switch(MO.ExtPulses) = MO.Amp(MO.ExtPulses);
%             MO.Offset(MO.ExtPulses) = MO.GetPhasePerc(Xcon);
%             % Set the torque to get turned off after the neuron fires if
%             % the off event is larger than 100% of osc. period
%             Overflow = MO.Offset(MO.ExtPulses)+MO.Duration(MO.ExtPulses)>1;
%             MO.Offset(MO.ExtPulses) = MO.Offset(MO.ExtPulses) - Overflow;
        end
        
        function PlotTorques(MO, tstep, mux)
            if nargin<3
                [Time,TorqueSig] = MO.GetTorqueSig();
            else
                [Time,TorqueSig] = MO.GetTorqueSig(tstep, mux);
            end

            plot(Time, TorqueSig);
        end
        
        function [Time, TorqueSig] = GetTorqueSig(MO, pstep, mux)
            % Check number of inputs to method
            if nargin<2
                pstep = 0.01;
            end
            if nargin<3
                mux = 1;
            end
            
            % Prepare empty variables
            Phase = (0:pstep:1)';
            Np = length(Phase);
            Time = Phase/MO.omega;
            if mux
                TorqueSig = zeros(Np,size(MO.OutM,1));
            else
                TorqueSig = zeros(Np,size(MO.OutM,2));
            end
            
            % Build torque signals
            for p = 1:MO.nPulses
                pStart = MO.Offset(p);
                pEnd = pStart + MO.Duration(p);
                if mux
                    j = find(MO.OutM(:,p)==1,1,'first');
                    TorqueSig(:,j) = TorqueSig(:,j) + ...
                        MO.Amp(p)*(Phase >= pStart & Phase < pEnd);
                else
                    TorqueSig(:,p) = ...
                        MO.Amp(p)*(Phase >= pStart & Phase < pEnd);
                end
            end
        end
        
        function diffs = CompareSigs(MO,NC2,pstep)
            if nargin<3
                pstep = 0.001;
            end
            [thisT,thisSig]=MO.GetTorqueSig(pstep,1);
            [~,otherSig]=NC2.GetTorqueSig(pstep,1);
            
            diffs = zeros(1,size(thisSig,2));
            try
                dSig = abs(thisSig-otherSig);
                for p = 1:size(thisSig,2)
                    diffs(p) = trapz(thisT,dSig(:,p));
                end
            catch err
                diffs = diffs + 1e3;
                disp(['Error comparing signals: ',err]);
            end
        end
            
    end
end