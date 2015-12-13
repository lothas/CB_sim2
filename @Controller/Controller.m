classdef Controller < handle & matlab.mixin.Copyable
    % Version 0.4 - 28/04/2014
    
    % 4 Leaky Integrate and Fire oscillators.
    % A pair of flexor and extensor LIF oscillator for each leg.
    
    properties
        % LIF parameters
        P_reset = 0;
        P_th = 1;
        
        P_LegE = 0.65; % timed leg extension

        % Oscillator's frequency
        omega = 1.106512566;
        omega0 = 1.106512566;
        
        % Initial phase for oscillator:
        P_0 = 0;
                
        stDim = 1; % state dimension
        nEvents = 2; % num. of simulation events
        
        % Controller Output
        nPulses = 0; % Overall number of pulses
        OutM;        % Output matrix (multiplies Switch vector)
        
        Amp = 0;     % Pulse amplitude in N
        Amp0 = 0;    % Base pulse amplitude in N
        Offset = 0;  % Defines beginning of pulse as % of neuron period
        Duration = 0;% Pulse duration as % of neuron period
        
        Switch;      % 0 when off, Amp when pulse is active
        ExtPulses;   % External pulses IDs
        
        % Feedback
        FBType = 2;  % 0 - no feedback
                     % 1 - single gain for omega and each joint
                     % 2 - individual gains for each pulse (not
                     % implemented?)
        lastPhi = 0;
        s_in = 0;    % Speed input - higher level input to control the
                     % desired walking speed
                     
        % Phase reset
        ExtP_reset = []; % set to a certain phase to use phase reset
                         % on foot contact
                         
        % Angular velocity impulses
        FBImpulse = 0; % 0 - no feedback
                       % 1 - add a constant value
                       % 2 - set ang. vel. to certain value
        AngVelImp = [];
        
        % Gains
        kOmega_u = 0.0;
        kOmega_d = 0.0;
        kTorques_u = 0;
        kTorques_d = 0;
        sOmega_f = 0.0;
        sOmega_s = 0.0;
        sTorques_f = 0;
        sTorques_s = 0;
        
        % Saturation
        MinSat; MaxSat;
        
        % Set keys
        SetKeys = {'P_reset','P_th','P_0','P_LegE','omega0',...
                   'nPulses','Amp0','Offset','Duration','AngVelImp',...
                   'FBType',...
                   'kOmega_u','kOmega_d','kTorques_u','kTorques_d',...
                   'sOmega_f','sOmega_s','sTorques_f','sTorques_s'};
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function [NC] = Controller(varargin)
            % Set torque parameters
            NC.nPulses = 4;
            NC.Amp0=[-13.6255 -4.6281 13.6255 0];
            NC.Amp=NC.Amp0;
            NC.OutM = [1 1 0 0; 0 0 1 1];
            
            NC.Offset=[0 0 0 0];
            
            NC.Duration=[0.1 0.4 0.1 0];
            
            NC.Switch=zeros(NC.nPulses,1);
            
            % Set adaptation parameters
            NC.kOmega_u = 0; %0.4579;
            NC.kOmega_d = 0; %0.9;
            NC.kTorques_u= 0; %[95 -443 95 0];
            NC.kTorques_d= 0; %[80 -295 80 0];
            
            NC.nEvents=1+NC.nPulses*2+1;
        end
        
        function [NC] = ClearTorques(NC)
            NC.nPulses = 0;
            NC.OutM = [0,0]';
            NC.Amp0 = [];
            NC.Amp = [];
            NC.Offset = [];
            NC.Duration = [];
            NC.Switch = 0;
            if NC.FBType == 2
                NC.kTorques_u = [];
                NC.kTorques_d = [];
            end
            NC.nEvents = 2;
            NC.ExtPulses = [];
        end
        
        function NC = Reset(NC,Phase)
            if nargin<2
                NC.Switch = 0*NC.Switch;
            else
                % Find which Torques should be active
                Start = NC.Offset;
                End = Start + NC.Duration;
                [Xperc] = NC.GetPhasePerc(Phase);
                On = Xperc >= Start & Xperc < End;
                NC.Switch = (NC.Amp.*On)';
            end 
        end
        
        function [Torques] = NeurOutput(NC)
            Torques = NC.OutM*NC.Switch;
        end

        function [per] = GetPeriod(NC)
            per = (NC.P_th-NC.P_reset)/NC.omega;
        end
        
        function [Xperc] = GetPhasePerc(NC,X)
            Xperc = (X-NC.P_reset)/(NC.P_th-NC.P_reset);
        end
        
        function diff = PhaseDiff(NC,ph1,ph2)
            diff = ph1 - ph2;
            % wrap it up
            Range = NC.P_th - NC.P_reset;
            diff(diff>Range/2) = diff(diff>Range/2) - Range;
            diff(diff<-Range/2) = diff(diff<-Range/2) + Range;
        end
        
        function [NC] = Adaptation(NC, Phi)
            if NC.FBType == 0
                % NO FEEDBACK
                NC.omega = NC.omega0;
                NC.Amp = NC.Amp0;
            else
                if abs(Phi-NC.lastPhi)>0.1
                    % Don't apply changes when the slope
                    % varies too abruptly
                    % (usually happens when the robot falls)
                    return
                end
                NC.lastPhi = Phi;
                
                NC.omega = NC.omega0 + ...
                    min(0,Phi)*NC.kOmega_d + ...    % Phi<0
                    max(0,Phi)*NC.kOmega_u;         % Phi>0
                NC.Amp = NC.Amp0 + ...
                    min(0,Phi)*NC.kTorques_d + ...  % Phi<0
                    max(0,Phi)*NC.kTorques_u;       % Phi>0
            end
            
            % Apply higher-level speed input
            NC.omega = NC.omega0 + ...
                min(0,NC.s_in)*NC.sOmega_s + ...    % s_in<0 (slow)
                max(0,NC.s_in)*NC.sOmega_f;         % s_in>0 (fast)
            NC.Amp = NC.Amp0 + ...
                min(0,NC.s_in)*NC.sTorques_s + ...  % s_in<0 (slow)
                max(0,NC.s_in)*NC.sTorques_f;       % s_in>0 (fast)
            
            % Apply saturation
            if ~isempty(NC.MinSat)
                NC.Amp = min(max(NC.Amp,NC.MinSat),NC.MaxSat);
            end
        end
        
        % %%%%%% % Derivative % %%%%%% %
        function [Xdot] = Derivative(NC, t, X) %#ok<INUSD>
            Xdot = NC.omega;
        end
        
        % %%%%%% % Events % %%%%%% %
        function [value, isterminal, direction] = Events(NC, X)
            value = ones(NC.nEvents,1);
            isterminal = ones(NC.nEvents,1);
            direction = -ones(NC.nEvents,1);
            
            % Check for firing neuron
            value(1) = NC.P_th - X;

            % Check for leg extension signal (for clearance)
%             value(2) = NC.P_LegE - X;
            
            % Check for switching on/off signal
            Xperc = NC.GetPhasePerc(X);
            value(3:2+NC.nPulses) = NC.Offset - Xperc;
            value(3+NC.nPulses:2+2*NC.nPulses) = ...
                NC.Offset + NC.Duration - Xperc;
        end
        
        function [NC,Xa] = HandleEvent(NC, EvID, Xb)
            Xa = Xb;
            switch EvID
                case 1
                    % Neuron fired
                    for i=1:NC.nPulses
                        if NC.Offset(i)==0;
                            % Turn on signal now
                            NC.Switch(i) = NC.Amp(i);
                        end
                    end
                    Xa = NC.P_reset; % reset phase
                case 2
                    % Extend the leg
                    % This is done from the simulation file
                case num2cell(3:2+NC.nPulses)
                    % Switch on signal
                    NC.Switch(EvID-2) = NC.Amp(EvID-2);
                case num2cell(3+NC.nPulses:2+2*NC.nPulses)
                    % Switch off signal
                    PulseID = EvID-(2+NC.nPulses);
                    NC.Switch(PulseID) = 0;
                    if any(PulseID == NC.ExtPulses)
                        % Set offset back to 200%
                        NC.Offset(PulseID) = 2;
                    end
            end
        end
        
        function [NC, Xmod, Xcon] = HandleExtFB(NC, Xmod, Xcon, Slope)
            % This function is called when the leg hits the ground
            if NC.FBType > 0
                % Perform adaptation based on terrain slope
                NC = NC.Adaptation(Slope);
            end
            
            if ~isempty(NC.ExtP_reset)
                % Perform a phase reset
                Xcon = NC.ExtP_reset;
                
                % Check if any event happens at ExtP_reset
                [value, it, dir] = NC.Events(Xcon); %#ok<NASGU,ASGLU>
                EvIDs = find(value == 0);
                for ev = 1:length(EvIDs)
                    [NC,Xcon] = NC.HandleEvent(EvIDs(ev),Xcon);
                end
            end
            
            switch NC.FBImpulse
                case 1 % 1 - add a constant value
                    Xmod(3:4) = Xmod(3:4) + NC.AngVelImp;
                case 2 % 2 - set ang. vel. to certain value
%                     delta = NC.AngVelImp - Xmod(3:4)
                    Xmod(3:4) = NC.AngVelImp;
            end
            
            % Activate external pulses
            NC.Switch(NC.ExtPulses) = NC.Amp(NC.ExtPulses);
            NC.Offset(NC.ExtPulses) = NC.GetPhasePerc(Xcon);
            % Set the torque to get turned off after the neuron fires if
            % the off event is larger than 100% of osc. period
            Overflow = NC.Offset(NC.ExtPulses)+NC.Duration(NC.ExtPulses)>1;
            NC.Offset(NC.ExtPulses) = NC.Offset(NC.ExtPulses) - Overflow;
        end
        
        function PlotTorques(NC,tstep,mux)
            if nargin<3
                [Time,TorqueSig]=NC.GetTorqueSig();
            else
                [Time,TorqueSig]=NC.GetTorqueSig(tstep,mux);
            end

            plot(Time,TorqueSig);
        end
        
        function [Time,TorqueSig]=GetTorqueSig(NC,pstep,mux)
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
            Time = Phase/NC.omega;
            if mux
                TorqueSig = zeros(Np,size(NC.OutM,1));
            else
                TorqueSig = zeros(Np,size(NC.OutM,2));
            end
            
            % Build torque signals
            for p = 1:NC.nPulses
                pStart = NC.Offset(p);
                pEnd = pStart + NC.Duration(p);
                if mux
                    j = find(NC.OutM(:,p)==1,1,'first');
                    TorqueSig(:,j) = TorqueSig(:,j) + ...
                        NC.Amp(p)*(Phase >= pStart & Phase < pEnd);
                else
                    TorqueSig(:,p) = ...
                        NC.Amp(p)*(Phase >= pStart & Phase < pEnd);
                end
            end
        end
        
        function diffs = CompareSigs(NC,NC2,pstep)
            if nargin<3
                pstep = 0.001;
            end
            [thisT,thisSig]=NC.GetTorqueSig(pstep,1);
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