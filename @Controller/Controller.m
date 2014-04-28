classdef Controller
    % Version 0.4 - 28/04/2014
    
    % 4 Leaky Integrate and Fire oscillators.
    % A pair of flexor and extensor LIF oscillator for each leg.
    
    properties
        % LIF parameters
        P_reset=0;
        P_th=1;
        
        P_LegE=0.65; % timed leg extension

        % Oscillator's frequency
        omega=1.106512566;
        omega0=1.106512566;
        
        % Initial phase for oscillator:
        P_0=0;
        
        NumEvents=0;
        
        % Controller Output
        NumTorques=4;   % Half applied at the ankle, half at the hip
        NumActJoints=2; % Active joints:
                        %    1 - Just hip actuation
                        %    2 - Ankle and Hip actuation
        
        NTorque=0;      % Output torque in N
        NTorque0=0;     % Base output torque in N
        NOffset=0;      % Defines beginning of pulse as % of neuron period
        NDuration=0;    % Pulse duration as % of neuron period
        
        NSwitch=0;
        tSon=0;
        tSoff=0;
        
        % Adaptation
        Adaptive=1;
        
        kOmega_up=0.4579;
        kOmega_down=0.9;
        kTorques_up=0;
        kTorques_down=0;
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function [NC] = Controller(varargin)
            % Set torque parameters
            NC.NTorque0=[-13.6255 -4.6281 13.6255 0];
            NC.NTorque=NC.NTorque0;
            
            NC.NOffset=[0 0 0 0];
            
            NC.NDuration=[0.1 0.4 0.1 0];
            
            NC.NSwitch=zeros(1,NC.NumTorques);
            NC.tSon=zeros(1,NC.NumTorques);
            NC.tSoff=zeros(1,NC.NumTorques);
            
            % Set adaptation parameters
            NC.kTorques_up=[95 -443 95 0];
            NC.kTorques_down=[80 -295 80 0];
            
            NC.NumEvents=1+NC.NumTorques*2+1;
        end
        
        function [Torques] = NeurOutput(NC)
            Torques=zeros(2,1); % [ Ankle; Hip ]
            
            switch NC.NumActJoints
                case 1
                    for i=1:NC.NumTorques
                        % Hip torque
                        Torques(2)=Torques(2)+NC.NSwitch(i)*NC.NTorque(i);
                    end
                case 2
                    half=NC.NumTorques/2;
                    for i=1:half
                        % Ankle torque
                        Torques(1)=Torques(1)+NC.NSwitch(i)*NC.NTorque(i);
                        % Hip torque
                        Torques(2)=Torques(2)+NC.NSwitch(half+i)*NC.NTorque(half+i);
                    end
            end
        end

        function [NC] = Adaptation(NC, X)
            if NC.Adaptive==0
                % NO FEEDBACK
                NC.omega=NC.omega0;
                NC.NTorque=NC.NTorque0;
            end
                
            if NC.Adaptive==1
                Phi=(X(1)+X(2))/2;

                if Phi<0
                    NC.omega=NC.omega0+NC.kOmega_down*Phi;
                    NC.NTorque=NC.NTorque0+NC.kTorques_down*Phi;
                else
                    NC.omega=NC.omega0+NC.kOmega_up*Phi;
                    NC.NTorque=NC.NTorque0+NC.kTorques_up*Phi;
                end
            end
        end
        
        % %%%%%% % Derivative % %%%%%% %
        function [Xdot] = Derivative(NC, t, X) %#ok<INUSD,MANUSL>
            Xdot=NC.omega;
        end
        
        % %%%%%% % Events % %%%%%% %
        function [value isterminal direction] = Events(NC, t, X)
            value=ones(NC.NumEvents,1);
            isterminal=ones(NC.NumEvents,1);
            direction=-ones(NC.NumEvents,1);
            
            % Check for firing neuron
            value(1)=NC.P_th-X(5);

            % Check for switching on signal
            value(2:1+NC.NumTorques)=NC.tSon-t;
            value(2+NC.NumTorques:1+2*NC.NumTorques)=NC.tSoff-t;
            
            % Check for leg extension signal (for clearance)
            value(end)=NC.P_LegE-X(5);
        end
        
        function [NC] = HandleEvents(NC, EvID, curTime)
            % EvID is a vector of length NumEvents, each item set to 1 if
            % event occurred (default is 0)
            if EvID(1)==1
                % Neuron fired
                % Set timers to turn on the torques
                NC.tSon=curTime+NC.NOffset/NC.omega;
                % Set timers to turn off the torques
                NC.tSoff=NC.tSon+NC.NDuration/NC.omega;
                
                for i=1:NC.NumTorques
                    if NC.NOffset(i)==0;
                        % Turn on signal now
                        NC.NSwitch(i)=1;
                    end
                end                        
            end
            for i=2:1+NC.NumTorques
                if EvID(i)==1
                    % Switch on signal
                    NC.NSwitch(i-1)=1;
                end
            end
            for i=2+NC.NumTorques:1+2*NC.NumTorques
                if EvID(i)==1
                    % Switch off signal
                    NC.NSwitch(i-(1+NC.NumTorques))=0;
                end
            end
            if EvID(end)==1
                % Extend the leg
                % This is done from the simulation file
            end
        end
        
        function [NC] = LoadParameters(NC,ControlParams)
            % Set CPG Parameters according to ControlParams input
            % ControlParams should include:
            % [ Omega,  leg extend phase, 
            %   Torque 1 strength, offset, duration,
            %   Torque 2 strength, offset, duration,
            %   ...
            %   Torque n strength, offset, duration]
            
            NC.NumTorques=(length(ControlParams)-2)/3;
            NC.NumEvents=1+NC.NumTorques*2+1;

            NC.omega0=ControlParams(1);
            NC.kOmega_up=0;
            NC.kOmega_down=0;

            NC.P_LegE=ControlParams(2);
            
            NC.NSwitch=zeros(1,NC.NumTorques);
            
            NC.NTorque0=zeros(1,NC.NumTorques);
            NC.NOffset=NC.NTorque0;
            NC.NDuration=NC.NTorque0;
            NC.kTorques_up=NC.NTorque0;
            NC.kTorques_down=NC.NTorque0;
            for i=1:NC.NumTorques
                NC.NTorque0(i)=ControlParams(3*i);
                NC.NTorque(i)=NC.NTorque0(i);
                NC.NOffset(i)=ControlParams(3*i+1);
                NC.NDuration(i)=ControlParams(3*i+2);
                NC.kTorques_up(i)=0;
                NC.kTorques_down(i)=0;
            end
        end
        
        function PlotTorques(NC)
            NumSteps=1000;
            TorqueSig=zeros(NC.NumTorques,NumSteps);

            Phase=linspace(0,1,NumSteps);
            for p=1:NumSteps
                for t=1:NC.NumTorques
                    if Phase(p)>=NC.NOffset(t) && Phase(p)<NC.NOffset(t)+NC.NDuration(t)
                        TorqueSig(t,p)=NC.NTorque(t);
                    end
                end
            end        

            plot(Phase,TorqueSig);
        end
        
        function [Time,TorqueSig]=GetTorqueSig(NC,tstep)
            Time=0:tstep:1/NC.omega0+tstep;
            NumSteps=length(Time);
            
            Phase=linspace(0,1,NumSteps);
            TorqueSig=zeros(NC.NumTorques,NumSteps);
            
            for p=1:NumSteps
                for t=1:NC.NumTorques
                    if Phase(p)>=NC.NOffset(t) && Phase(p)<NC.NOffset(t)+NC.NDuration(t)
                        TorqueSig(t,p)=NC.NTorque(t);
                    end
                end
            end        
        end
    
    end
end