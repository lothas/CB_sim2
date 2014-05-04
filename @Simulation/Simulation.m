classdef Simulation < handle
    % Version 0.1 - 04/05/2014
    % This simulation integrates a system over time until
    % an event occurs, then it performs some calculations
    % and continues integrating until the next event or
    % it runs out of time
    
    properties
        Mod; % Model
        Con; % Controller
        Env; % Environment
    
        % State params
        stDim; ModCo; ConCo;
        % Event params
        nEvents; ModEv; ConEv;
    
        % Simulation parameters
        IC;
        infTime;
        tstart; tstep; tend; tspan;
        
        % Performance tracking
        CurSpeed; StepsTaken;
                
        % Rendering params
        Graphics = 1;
        Fig = 0; Once; StopSim;
        FigWidth; FigHeight; AR;
        % Environment display
        FlMin; FlMax; HeightMin; HeightMax;
        
        Follow = 1; % Follow model with camera
        % COM transformation
        tCOM; COMx0; COMy0;
        % Time display
        hTime; TimeStr = ['t = %.2f s\nOsc.=%.3f\n',...
                         'Slope = %.2f ',char(176)','\nSpeed = %s'];
        % Torques display
        nOuts; nTsteps;
        Ttime; Thold; Tbase; Tscale;
        hTorques;
        Colors = {[1 0 0],[0 0 1],[0 1 0],[0 0 0]};
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function sim = Simulation(varargin)
            switch nargin
                case 3
                    sim.Mod = varargin{1};
                    sim.Con = varargin{2};
                    sim.Env = varargin{3};
                otherwise
                    sim.Mod = CompassBiped();
                    sim.Con = Controller();
                    sim.Env = Terrain();
            end            
        end
        
        function sim = SetTime(sim,tstart,tstep,tend)
            if nargin~=4
                error(['Set time expects 3 input arguments',...
                    ' but was provided with ',num2str(nargin)]);
            end
            sim.tstart = tstart;
            sim.tstep = tstep;
            if isnumeric(tend)
                if tend<=tstart+tstep
                    error('tend is too close to tstart');
                else
                    sim.tend = tend;
                end
                sim.infTime = 0;
            else
                if strcmp(tend,'inf')
                    % Simulation will run for indefinite time
                    sim.infTime = 1;
                    sim.tend = 10;
                end
            end
        end
        
        function [Xt] = Derivative(sim,t,X)
            sim.Mod.Torques=sim.Con.NeurOutput();

            Xt = [sim.Mod.Derivative(t,X(sim.ModCo));
                  sim.Con.Derivative(t,X(sim.ConCo))];
        end

        function [value, isterminal, direction] = Events(sim, t, X) %#ok<INUSL>
            value = zeros(sim.nEvents,1);
            isterminal = ones(sim.nEvents,1);
            direction = zeros(sim.nEvents,1);

            % Call model event function
            [value(sim.ModEv), isterminal(sim.ModEv), direction(sim.ModEv)] = ...
                sim.Mod.Events(X(sim.ModCo), sim.Env);
            % Call controller event function
            [value(sim.ConEv), isterminal(sim.ConEv), direction(sim.ConEv)] = ...
                sim.Con.Events(X(sim.ConCo));
        end

        function StopButtonCb(sim, hObject, eventdata, handles) %#ok<INUSD>
            if sim.StopSim == 0
                sim.StopSim=1;
                set(hObject,'String','Close Window');
            else
                close(sim.Fig)
            end
        end  % StopButtonCallback
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ControlParams=[1.266646286756660, 0.597267606409947,...
%                -7.384173104761088, 0.126800300703305, 0.072267297914989,...
%                 5.191278913625541, 0.166498744328709, 0.053696788954836];
%                     
% SimType=1;
% 
% if SimType==1
%     % Terrain parameters
%     FloorType=0;
%     InitSlope=0;
%     EndSlope=0;
%     Adaptive=0;
% else
%     % Terrain parameters
%     FloorType=2;
%     InitSlope=0;
%     Adaptive=1;
%     
%     % Start the slope after giving the robot
%     % 10 seconds to reach steady state
%     Floor.start_x=AvgVel*10;
%     
%     MinSlope=0;
%     MaxSlope=0;
% end
% 
%     function Xout=ImpactHandle(X)
%         SS=1;
% 
%         % Update min/max slope
%         if SimType==2
%             % Update only if the robot is really walking
%             % not if it "leaped" while falling
%             if abs(X(2)-X(1))<0.4*pi
%                 % Get slope on spot of "hind" leg
%                 Slope=Floor.SurfSlope(min(Robot.xS,xNS))*180/pi;
%                 if Slope>MaxSlope
%                     MaxSlope=Slope;
%                 end
%                 if Slope<MinSlope
%                     MinSlope=Slope;
%                 end
%             end
%         end
%     end
% 
%     function StopButtonCallback(t,X) %#ok<INUSD>
%         StopSim=1;
%         close all
%         dispLength='%5.7f';
%         disp(['InitCond=[',num2str(IC_Store(1,1),dispLength),', ',num2str(IC_Store(2,1),dispLength),', ',...
%                         num2str(IC_Store(3,1),dispLength),', ',num2str(IC_Store(4,1),dispLength),', ',...
%                         num2str(IC_Store(5,1),dispLength),'];',10]);
%     end  % StopButtonCallback

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Output the end condition reached
% EndReached=0;
% EndText='Reached end of tspan';
% 
% % Variables for End Cond 1: Convergence
% Num_IC=10;
% IC_Store=zeros(5,Num_IC);
% 
% 
% %% Simulation
% 
% options=odeset('MaxStep',tstep/10,'RelTol',.5e-7,'AbsTol',.5e-8, 'OutputFcn', @RealTimePlot, 'Events', @EventFunc);
% [TTemp,XTemp,TE,YE,IE]=ode45(@Derivative,tspan,InitCond,options); %#ok<ASGLU>
% 
% T=TTemp;
% X=[ones(length(TTemp),1)*Robot.xS, ones(length(TTemp),1)*Robot.yS, XTemp];
% % Save torques
% TorquesTemp=zeros(2,length(TTemp));
% for j=1:length(TTemp)
%     TorquesTemp(:,j)=CPG.NeurOutput();
% end
% TorquesP=TorquesTemp;
% 
% % Simulation information
% StepsTaken=0;
% 
% while TTemp(end)<tspan(end-1) && StopSim==0
%     SS=0;
%     Xnew=XTemp(end,1:4);
%     OscReset=0;
%     
%     CPGEvents=zeros(1,CPG.NumEvents);
%     
%     for i=1:size(IE,1)
%         % Check Compass Biped events
%         switch IE(i)
%             case 1 % Foot contact
%                 Xnew=ImpactHandle(XTemp(end,:));
%                 
%                 % Check "End-Game" conditions
%                 if abs(Xnew(2)-Xnew(1))<0.00001
%                     EndReached=2;
%                     EndText='Step length too small';
%                     StopSim=1;
%                     break;
%                 end
%                 
%                 if Robot.LegShift>0
%                     % Robot hit the ground before extending the leg
%                     EndReached=3;
%                     EndText='Robot hit the ground before extending the leg';
%                     StopSim=1;
%                     break;
%                 end
%             case 2 % Hip height too low
%                 EndReached=1;
%                 EndText='Hip height too low';
%                 StopSim=1;
%                 break;
%             otherwise
%                 % Add event to CPG events
%                 CPGEvents(IE(i)-Robot.NumEvents)=1;
% 
%                 if IE(i)==Robot.NumEvents+1
%                     % Oscillator fired
%                     % Part of this event is handled here
%                     OscReset=CPG.P_reset-XTemp(end,5);
%                     Robot.LegShift=Robot.Clearance;
%                 end
% 
%                 if IE(i)==TotEvents
%                     % Leg extension timer
%                     % This event is handled here
%                     Robot.LegShift=0;
%                 end
%         end
%         CPG=CPG.HandleEvents(CPGEvents, TTemp(end));
%     end
%         
%     % Check hip height
%     [HipPosx,HipPosy]=Robot.GetPos(Xnew,'Hip'); %#ok<ASGLU>
%     if HipPosy-0.7*Robot.L<0            
%         EndReached=1;
%         EndText='Hip height too low';
%         StopSim=1;
%         break;
%     end
%     
%     IC=[Xnew, XTemp(end,5)+OscReset];
%         
%     if SS==1
%         IC_Store(:,2:end)=IC_Store(:,1:end-1);
%         IC_Store(:,1)=IC';
%     end
%         
%     % Continue simulation
%     if Graphics==2 || SimType==2
%         % Extend simulation
%         tend=tend+TTemp(end)-TTemp(1);
%     end
%     tspan=TTemp(end):tstep:tend;
%     [TTemp,XTemp,TE,YE,IE]=ode45(@Derivative,tspan,IC,options); %#ok<ASGLU>
%     
%     T=[T; TTemp]; %#ok<AGROW>
% %     if Robot.Support==Robot.Left
% %         XTemp2=[XTemp(:,2), XTemp(:,1), XTemp(:,4), XTemp(:,3), XTemp(:,5:end)];
% %         X=[X; ones(length(TTemp),1)*Robot.xS, ones(length(TTemp),1)*Robot.yS, XTemp2]; %#ok<AGROW>
% %     else        
% %         X=[X; ones(length(TTemp),1)*Robot.xS, ones(length(TTemp),1)*Robot.yS, XTemp]; %#ok<AGROW>
% %     end
%     X=[X; ones(length(TTemp),1)*Robot.xS, ones(length(TTemp),1)*Robot.yS, XTemp]; %#ok<AGROW>
%     
%     % Save torques
%     TorquesTemp=zeros(2,length(TTemp));
%     for j=1:length(TTemp)
%         TorquesTemp(:,j)=CPG.NeurOutput();
%     end
%     TorquesP=[TorquesP, TorquesTemp]; %#ok<AGROW>
% end
% 
% if TTemp(end)>=tspan(end-1)
%     EndReached=3;
%     EndText='Reached end of time span';
% end
% 
% switch nargout
%     case 0
%     	varargout={};
%         if Graphics~=2
%             disp(EndText);
%             dispLength='%5.7f';
%             disp(['InitCond=[',num2str(IC_Store(1,1),dispLength),', ',num2str(IC_Store(2,1),dispLength),', ',...
%                             num2str(IC_Store(3,1),dispLength),', ',num2str(IC_Store(4,1),dispLength),', ',...
%                             num2str(IC_Store(5,1),dispLength),'];',10]);
%         end
%     case 1
%         varargout={EndReached};
%         disp(EndText);
%         dispLength='%5.7f';
%         disp(['InitCond=[',num2str(IC_Store(1,1),dispLength),', ',num2str(IC_Store(2,1),dispLength),', ',...
%                         num2str(IC_Store(3,1),dispLength),', ',num2str(IC_Store(4,1),dispLength),', ',...
%                         num2str(IC_Store(5,1),dispLength),'];',10]);
%     case 2
%         if SimType==1
%             varargout={EndReached,IC_Store(:,1)};
%         else
%             varargout={EndReached,max(abs(MinSlope),abs(MaxSlope))};
%         end
%     case 3
%         varargout={IC_Store,1,1};
%     case 4
%         varargout={EndReached, EndText, 1, IC_Store(:,1)};
%     case 5
%         varargout={EndReached, EndText, T, X, TorquesP};
%     otherwise
%         varargout={};
% end
% end

