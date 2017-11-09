classdef Matsuoka < handle & matlab.mixin.Copyable
    % Version 0.1 - 04/01/2016
    
%%    % Matsuoka oscillator class
%     
% *) 'y'- is the output from all of the neurons

% *) 'OutM' - is the matrix which convert the CPG output ('y) to joints
%   torques. for example, is we have 4 neurons:
%   
%     joints_Torques = OutM*y = |1,-1,0, 0|  * |y1| = |y1 - y2| = |tau_ankle|
%                               |0, 0,1,-1|    |y2|   |y3 - y4|   |tau_hip  |
%                                              |y3|
%                                              |y4|
%     
%     Note: the ankle torque always comes first.
% 
% *) every step we also need to switch the sign of the CPG's hip torque.
%       #) In the simulation the stance and swing leg are always switching
%           at ground impact.
%       #) The hip torque is defined in the simulation as the torque
%           between the stance leg and the swing leg.
%       #) we want the Matsuoka CPG to represent the actual torque of a
%           real robot actuator. i.e. the CPG output need to represent the
%           torque between the left leg to the right, regardless of which one
%           of them is 'stance' or 'swing'.
%     So, every time that in the model the stance and the swing are
%     changing, we also chage the sign of the torque by changing the last
%     row of 'OutM' which represents always the hip torque.
%   Every ground impact we also multiply 'OutM' by [1,0;0,-1];
% 
% *) Important: if we use two ankle joints, then we need to switch between
%           which pair of neurons go to the stance leg every step.
%           'OutM' is defined:
%           |1,-1,0, 0,0, 0| = every even step
%           |0, 0,0, 0,1,-1|
%           |0, 0,1,-1,0, 0| = every odd step
%           |0, 0,0, 0,1,-1|
% 


%%
    properties
        name = 'Matsuoka'
        
%         Parameters for 4 nueons CPG
        tau0 = 1; tau = 1;
        tav0 = 1; tav = 1;
        tau_ratio = 12;  %5;
        beta = 0.1;
        u0 = 1;
        win = [];%[0, 3; 1, 0];
        wex = [];
        W = [];
     
        stDim = 4; % state dimension
        nEvents = 1; % num. of simulation events
        
        % Controller Output
        startup_t = 0;
        nPulses = 1;
        OutM = [0, 0; 1 -0.1];
        Amp0 = 1.5;%10
        Amp = 1.5;%10
        ExtPulses = [];
        
        % For case with two distinct ankle joints
        % every step a different pair of CPG neurons is chosen to actuate
        % the ankle.
        jointSelM = 1; %diag([1,1,0,0,1,1]);
        
        % if 'ture' then use 6neurons and actuate different ankles
        %   if 'false then only actuate 'hip'
        twoAnkleFlag = 0;
        
        
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
            'tau', 'tau_u', 'tav', 'tau_v', 'tau_r', '\tau_r', ...
            'tau_ratio', '\tau_ratio', ...
            'beta', 'win', 'wex', 'weights', 'amp0', 'amp', ...
            'fbtype', 'feedback', 'fb', ...
            'ks_tau', 'speed_tau', 'tau_speed_gain', 'ks_\tau', ...
            'ks_out', 'speed_out', 'torque_speed_gain', 'ks_c',...
            '2neuron_symm_weights','2neuron_general_weights',...
            'amp_2n_dif_inputs','amp_2n_same_inputs','ks_c_2n_symm','ks_c_2n_general',...
            '4neuron_symm_weights','amp_4n_symm','4neuron_taga_like','ks_c_4n_symm',...
            'amp_same4each_joint','ks_c_same4each_joint',...
            'amp_6n_symm','6neuron_taga_like'};
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function [MO] = Matsuoka(varargin)
            
%             MO.nEvents = ;
        end
        
        function MO = Reset(MO,Phase)
            
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
        
        function MO = SetAnkles_selection(MO, ankleNum)
            % If we use two ankles the set 'jointSelM' that every step is
            % chaging via event to be either:
            %   diag([1,1,0,0,1,1]) - when the left leg is stance
            % or,
            %   diag([0,0,1,1,1,1]) - when the right leg is stance
            if ankleNum == 2
                MO.jointSelM = diag([1,1,0,0,1,1]);
                MO.twoAnkleFlag = 1;
            elseif ankleNum == 1
                MO.jointSelM = 1;
                MO.twoAnkleFlag = 0;
            end
        end 
        
        function [Torques] = Output(MO, t, MOX, ~)
            y = max(MOX(1:2*MO.nPulses,:),0);
%             y = diag(MO.Amp)*max(MOX(1:2*MO.nPulses,:),0);
            try
%                 Torques = MO.OutM*y;
                Torques = MO.OutM*MO.jointSelM*y;
            catch
                disp('unable to calculate torques')
            end

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
            MO.tau = MO.tau0 + MO.ks_tau*MO.s_in;
            MO.tau = max(MO.tau, 0.01); % Tau has to be > 0
            MO.tav = MO.tau_ratio*MO.tau;
            MO.Amp = MO.Amp0 + MO.ks_out*MO.s_in;
            MO.Amp = max(MO.Amp, 0.01); % Amp has to be > 0
            % Update the weight matrix
            MO.W = (diag(1./MO.Amp)*(diag(MO.Amp)*(MO.win+MO.wex))')';
%             MO.W = min(max(MO.W, -10),100); % Bound the weights
        end
        
        % %%%%%% % Derivative % %%%%%% %
        function [Xdot] = Derivative(MO, ~, X)
            % X = [u_i; v_i];
            u = X(1:2*MO.nPulses,:);
            v = X(2*MO.nPulses+1:end,:);
            y = max(u,0);
            
            udot = 1/MO.tau*(-u - MO.beta*v + MO.Amp - MO.W*y);
            vdot = 1/MO.tav*(-v+y);
            
            Xdot = [udot;
                    vdot];
        end
        
        % %%%%%% % Events % %%%%%% %
        function [value, isterminal, direction] = Events(MO, X)
            value = ones(MO.nEvents,1);
            isterminal = ones(MO.nEvents,1);
            direction = -ones(MO.nEvents,1);

        end
        
        function [MO,Xa] = HandleEvent(MO, EvID, Xb)
            Xa = Xb;
        end
        
        function [MO, Xmod, Xcon] = HandleExtFB(MO, Xmod, Xcon, Slope)
            % This function is called when the leg hits the ground
%             if MO.FBType > 0
%                 % Perform adaptation based on terrain slope
                MO = MO.Adaptation(Slope);
%             end
                
                if MO.twoAnkleFlag
                    % define 3 pairs of neurons, thus both ankles have a
                    %   different signal.
                    MO.jointSelM = abs(MO.jointSelM - diag([1,1,1,1,0,0]));
                end
                
                % % %reverse the hip torque after each ground impact
                MO.OutM = [1,0; 0, -1]*MO.OutM;
                
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