classdef Simulation
    % Version 0.1 - 27/04/2014
    % This simulation integrates a system over time until
    % an event occurs, then it performs some calculations
    % and continues integrating until the next event or
    % it runs out of time
    
    properties
        Mod; % Model
        Con; % Controller
        Env; % Environment
        
        % Time
        tstart = 0; tstep = 0.05; tend = 30; tspan;
        
        stDim = 0;
        ModSt = 0; ConSt = 0;
        nEvents = 0;
        ModEv = 0; ConEv = 0;
        
        % External functions handles
        InitFcn = 0;
        StepFcn = 0;
        EventFcn = 0;
        EndFcn = 0;
        
        % Rendering params
        Graphics = 1;
        StopSim = 0;
        Follow = 1; % Follow model with camera
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function sim = Simulation(varargin)
            switch nargin
                case 3
                    sim.Mod = varargin{1};
                    sim.Con = varargin{2};
                    sim.Env = varargin{3};
            end            
        end
        
        function sim = Init(sim, MP, CP, EP)
            % Receives model parameters (MP), controller
            % parameters (CP) and environment parameters (EP)
            % and initializes all the objects
            
            % Initialize robot
            sim.Mod = CompassBiped(MP);
            
            % Initialize the controller (CPG)
            sim.Con = Controller();
            sim.Con=sim.Con.LoadParameters(CP);
            
            % Adapt CPG (if adaptive)
            sim.Con.Adaptive=Adaptive;
            sim.Con=sim.Con.Adaptation([0 2*Floor.SurfSlope(Robot.xS) 0 0 0]);
            
            % Initialize the terrain
            sim.Env = Terrain(EP);
            
            % Set simulation time span
            sim.tspan=sim.tstart:sim.tstep:sim.tend;
            
            % Set the state IDs
            ModDim = sim.Mod.stDim;
            sim.ModSt = 1:ModDim;
            ConDim = sim.Con.stDim;
            sim.ConSt = ModDim+1:ModDim+ConDim;
            sim.stDim = ModDim+ConDim;
            
            % Set total number of simulation events
            nMEv = sim.Mod.nEvents;
            sim.ModEv = 1:nMEv;
            nCEv = sim.Con.nEvents;
            sim.ConEv = nMEv+1:nMEv+nCEv;
            sim.nEvents=nMEv+nCEv;
    
            if isa(sim.InitFcn,'function_handle')
                % Call external function
                sim = sim.InitFcn(sim);
            end
        end
        
        function sim = Step(sim, t, X)
            
            if isa(sim.StepFcn,'function_handle')
                % Call external function
                sim = sim.StepFcn(sim);
            end
        end
        
        function sim = Events(sim, evID, t, X)
            if isa(sim.EventFcn,'function_handle')
                % Call external function
                sim = sim.EventFcn(sim);
            end
        end
        
        function sim = End(sim)
            if isa(sim.EndFcn,'function_handle')
                % Call external function
                sim = sim.EndFcn(sim);
            end
        end
    end
end

function [ varargout ] = Simulation( varargin )
% Version 0.7 - 18/12/2012

SimType=1;

switch nargin
    case 4
        ControlParams=varargin{1};
        NumActJoints=varargin{2};
        tend=varargin{3};
        Graphics=varargin{4};
    case 5
        ControlParams=varargin{1};
        NumActJoints=varargin{2};
        SlopeCurve=varargin{3};
        AvgVel=varargin{4};
        Graphics=varargin{5};
        SimType=2;
        
        tend=5; % Time will be extended automatically
    otherwise
        clc
        ControlParams=[1.266646286756660, 0.597267606409947,...
                       -7.384173104761088, 0.126800300703305, 0.072267297914989,...
                        5.191278913625541, 0.166498744328709, 0.053696788954836];
       NumActJoints=2; 
        tend=30;
end

% Start the robot from and idle position
InitCond=[0 0 0 0 0];

% Initialize robot and shorten the swing leg for clearance
Robot = CompassBiped();
Robot.LegShift=Robot.Clearance;

CurSpeed='Computing...';

% Initialize the controller (CPG)
CPG = Controller();
CPG.NumActJoints=NumActJoints;

CPG=CPG.LoadParameters(ControlParams);
    
% Set total number of simulation events
TotEvents=Robot.NumEvents+CPG.NumEvents;

if SimType==1
    % Terrain parameters
    FloorType=0;
    InitSlope=0;
    EndSlope=0;
    Adaptive=0;

    % Initialize the terrain
    Floor = Terrain(FloorType,InitSlope,EndSlope);
else
    % Terrain parameters
    FloorType=2;
    InitSlope=0;
    Adaptive=1;
    
    % Initialize the terrain
    Floor = Terrain(FloorType,InitSlope,SlopeCurve);
    
    % Start the slope after giving the robot
    % 10 seconds to reach steady state
    Floor.start_x=AvgVel*10;
    
    MinSlope=0;
    MaxSlope=0;
end

% Adapt CPG (if adaptive)
CPG.Adaptive=Adaptive;
CPG=CPG.Adaptation([0 2*InitSlope*pi/180 0 0 0]);

% Set simulation time span
tstart=0;
tstep=0.05;
tspan=tstart:tstep:tend;

% Call oscillator fire event
Events=zeros(CPG.NumEvents);
Events(1)=1;
CPG=CPG.HandleEvents(Events, tstart);

% Set display parameters
global StopSim;
StopSim=0;

% Scene parameters
if nargin<3
    Graphics=1;
end
Follow=1;

L=Robot.L;
if Follow
    [COMx,COMy]=Robot.GetPos(InitCond,'COM');
    [COMy,A]=Floor.Surf(COMx);
    FlMin=COMx-2*L;
    FlMax=COMx+2*L;
    HeightMin=COMy-1.5*L;
    HeightMax=COMy+1.5*L;
else
    FlMin=-2*L;
    FlMax=2*L;
    HeightMin=-1.5*L;
    HeightMax=1.5*L;
end

scrsz=0;
FigWin=0;
tCOM=0;
COMx0=0;
COMy0=0;
TimeDisp=0;
hTankle=0;
hThip=0;
Once=1;

    function [Xdot] = Derivative(t, X)
        Robot.Torques=CPG.NeurOutput();
            
        Xdot=[ Robot.Derivative(t, X);
        	    CPG.Derivative(t, X) ];
    end
    
    function [value isterminal direction] = EventFunc(t, X)
        value=ones(TotEvents,1);
        isterminal=ones(TotEvents,1);
        direction=-ones(TotEvents,1);
        
        % Check for foot contact / detachment
        % Check for swing stop
        value(1:Robot.NumEvents)=Robot.Events(X, Floor);
        
        % Check for firing neuron
        % Check for switching signal
        value(Robot.NumEvents+1:TotEvents)=CPG.Events(t, X);
        
        if StopSim==1
            value(:) = 0;
            direction(:) = 0;
        end
    end

    function Xout=ImpactHandle(X)
        [xNS, yNS]=Robot.GetPos(X,'NS');
        
        % Calculate Impact
        Xout=Robot.CalcImpact(X);

        % Update support foot
        if Robot.Support==Robot.Left
            Robot.Support=Robot.Right;
        else
            if Robot.Support==Robot.Right
            Robot.Support=Robot.Left;
            end
        end
        SS=1;

        % Update min/max slope
        if SimType==2
            % Update only if the robot is really walking
            % not if it "leaped" while falling
            if abs(X(2)-X(1))<0.4*pi
                % Get slope on spot of "hind" leg
                Slope=Floor.SurfSlope(min(Robot.xS,xNS))*180/pi;
                if Slope>MaxSlope
                    MaxSlope=Slope;
                end
                if Slope<MinSlope
                    MinSlope=Slope;
                end
            end
        end

        % Update support foot position
        Robot.yS=yNS;
        Robot.xS=xNS;

        % Adapt CPG (if adaptive)
        CPG=CPG.Adaptation(X);

        % Count step
        StepsTaken=StepsTaken+1;

        % Get Speed
        CurSpeed=Robot.GetStepLength(X)*CPG.omega;   
    end

    function status = RealTimePlot(t,X,flag)
        if strcmp(flag,'done')
            % Finish simulation
            return
        end
               
        if Graphics && StopSim==0
            if strcmp(flag,'init')
                if Once
                    % Initialize
                    if Graphics==1
                        scrsz = get(0, 'ScreenSize');
                        FigWin=figure();
                        set(FigWin,'Position', [150 100 scrsz(3)-500 scrsz(4)-300]);
                    end
                    set(gca,'LooseInset',get(gca,'TightInset')*2)

                    hold on
                    axis([FlMin FlMax HeightMin HeightMax]);
                    axis equal

                    [COMx0,COMy0]=Robot.GetPos(X,'COM'); 
                    [COMy0,A]=Floor.Surf(COMx); %#ok<SETNU>
                    tCOM = hgtransform('Parent',gca);

                    % Draw Robot
                    Robot=Robot.Render(X);
                    % Draw Floor
                    Floor=Floor.Render(FlMin,FlMax);

                    % Display time
                    TimeDisp=text(-1.2*L, 1.4*L, sprintf('t=%.2f',t),'HorizontalAlignment','left','Parent',tCOM);

                    % Display Torques
                    hTankle=text(FlMax-0.65*L, 0.13*L, '', 'BackgroundColor', [1 0.7 0.7],'Parent',tCOM);
                    hThip=text(FlMax-0.65*L, 0.26*L, '', 'BackgroundColor', [1 1 0.7],'Parent',tCOM);

                    if Graphics==1
                        uicontrol('Style', 'pushbutton', 'String', 'Stop Simulation',...
                            'Position', [scrsz(3)-700 scrsz(4)-360 100 30],...
                            'Callback', @StopButtonCallback);
                    end
                    if Graphics==2
                        uicontrol('Style', 'pushbutton', 'String', 'Stop Simulation',...
                            'Units','Normalized','Position', [0.9 0.89 0.08 0.05],...
                            'Callback', @StopButtonCallback2);
                    end

                    Once=0;
                end
                
                return
            end

            if Follow
                [COMx,COMy]=Robot.GetPos(X,'COM');
                [COMy,A]=Floor.Surf(COMx);
                FlMin=COMx-1.25*L;
                FlMax=COMx+1.25*L;
                HeightMin=COMy-0.5*L;
                HeightMax=COMy+1.5*L;
                
                TCOMx = makehgtform('translate',[COMx-COMx0 COMy-COMy0 0]);
                set(tCOM,'Matrix',TCOMx);
                
                axis([FlMin FlMax HeightMin HeightMax]);
            end
                
            % Draw Robot
            Robot.Render(X);
            % Draw Floor
            Floor=Floor.Render(FlMin,FlMax);
   
            % Display time
            if ischar(CurSpeed)
                set(TimeDisp,'String',sprintf(['t=%.2f - Osc.=%.3f\n',...
                                               'Slope=%.2f\n',...
                                               'Speed= %s'],t(1),X(5),Floor.SurfSlope(Robot.xS)*180/pi,CurSpeed));
            else
                set(TimeDisp,'String',sprintf(['t=%.2f - Osc.=%.3f\n',...
                                               'Slope=%.2f\n',...
                                               'Speed= %.3f m/s (%.3f km/h)'],t(1),X(5),Floor.SurfSlope(Robot.xS)*180/pi,CurSpeed,CurSpeed*3.6));
            end
            
            % Display Torques
            Torques=CPG.NeurOutput();
            if Torques(1)~=0
                set(hTankle,'String', 'Flexor');
            else
                set(hTankle,'String', '');
            end
            if Torques(2)~=0
                set(hThip,'String', 'Extensor');
            else
                set(hThip,'String', '');
            end

            drawnow
        end
        
        status=StopSim;
    end

    function StopButtonCallback(t,X) %#ok<INUSD>
        StopSim=1;
        close all
        dispLength='%5.7f';
        disp(['InitCond=[',num2str(IC_Store(1,1),dispLength),', ',num2str(IC_Store(2,1),dispLength),', ',...
                        num2str(IC_Store(3,1),dispLength),', ',num2str(IC_Store(4,1),dispLength),', ',...
                        num2str(IC_Store(5,1),dispLength),'];',10]);
    end  % StopButtonCallback

    function StopButtonCallback2(t,X) %#ok<INUSD>
        StopSim=1;
    end  % StopButtonCallback

% Output the end condition reached
EndReached=0;
EndText='Reached end of tspan';

% Variables for End Cond 1: Convergence
Num_IC=10;
IC_Store=zeros(5,Num_IC);


%% Simulation

options=odeset('MaxStep',tstep/10,'RelTol',.5e-7,'AbsTol',.5e-8, 'OutputFcn', @RealTimePlot, 'Events', @EventFunc);
[TTemp,XTemp,TE,YE,IE]=ode45(@Derivative,tspan,InitCond,options); %#ok<ASGLU>

T=TTemp;
X=[ones(length(TTemp),1)*Robot.xS, ones(length(TTemp),1)*Robot.yS, XTemp];
% Save torques
TorquesTemp=zeros(2,length(TTemp));
for j=1:length(TTemp)
    TorquesTemp(:,j)=CPG.NeurOutput();
end
TorquesP=TorquesTemp;

% Simulation information
StepsTaken=0;

while TTemp(end)<tspan(end-1) && StopSim==0
    SS=0;
    Xnew=XTemp(end,1:4);
    OscReset=0;
    
    CPGEvents=zeros(1,CPG.NumEvents);
    
    for i=1:size(IE,1)
        % Check Compass Biped events
        switch IE(i)
            case 1 % Foot contact
                Xnew=ImpactHandle(XTemp(end,:));
                
                % Check "End-Game" conditions
                if abs(Xnew(2)-Xnew(1))<0.00001
                    EndReached=2;
                    EndText='Step length too small';
                    StopSim=1;
                    break;
                end
                
                if Robot.LegShift>0
                    % Robot hit the ground before extending the leg
                    EndReached=3;
                    EndText='Robot hit the ground before extending the leg';
                    StopSim=1;
                    break;
                end
            case 2 % Hip height too low
                EndReached=1;
                EndText='Hip height too low';
                StopSim=1;
                break;
            otherwise
                % Add event to CPG events
                CPGEvents(IE(i)-Robot.NumEvents)=1;

                if IE(i)==Robot.NumEvents+1
                    % Oscillator fired
                    % Part of this event is handled here
                    OscReset=CPG.P_reset-XTemp(end,5);
                    Robot.LegShift=Robot.Clearance;
                end

                if IE(i)==TotEvents
                    % Leg extension timer
                    % This event is handled here
                    Robot.LegShift=0;
                end
        end
        CPG=CPG.HandleEvents(CPGEvents, TTemp(end));
    end
        
    % Check hip height
    [HipPosx,HipPosy]=Robot.GetPos(Xnew,'Hip'); %#ok<ASGLU>
    if HipPosy-0.7*Robot.L<0            
        EndReached=1;
        EndText='Hip height too low';
        StopSim=1;
        break;
    end
    
    IC=[Xnew, XTemp(end,5)+OscReset];
        
    if SS==1
        IC_Store(:,2:end)=IC_Store(:,1:end-1);
        IC_Store(:,1)=IC';
    end
        
    % Continue simulation
    if Graphics==2 || SimType==2
        % Extend simulation
        tend=tend+TTemp(end)-TTemp(1);
    end
    tspan=TTemp(end):tstep:tend;
    [TTemp,XTemp,TE,YE,IE]=ode45(@Derivative,tspan,IC,options); %#ok<ASGLU>
    
    T=[T; TTemp]; %#ok<AGROW>
%     if Robot.Support==Robot.Left
%         XTemp2=[XTemp(:,2), XTemp(:,1), XTemp(:,4), XTemp(:,3), XTemp(:,5:end)];
%         X=[X; ones(length(TTemp),1)*Robot.xS, ones(length(TTemp),1)*Robot.yS, XTemp2]; %#ok<AGROW>
%     else        
%         X=[X; ones(length(TTemp),1)*Robot.xS, ones(length(TTemp),1)*Robot.yS, XTemp]; %#ok<AGROW>
%     end
    X=[X; ones(length(TTemp),1)*Robot.xS, ones(length(TTemp),1)*Robot.yS, XTemp]; %#ok<AGROW>
    
    % Save torques
    TorquesTemp=zeros(2,length(TTemp));
    for j=1:length(TTemp)
        TorquesTemp(:,j)=CPG.NeurOutput();
    end
    TorquesP=[TorquesP, TorquesTemp]; %#ok<AGROW>
end

if TTemp(end)>=tspan(end-1)
    EndReached=3;
    EndText='Reached end of time span';
end

switch nargout
    case 0
    	varargout={};
        if Graphics~=2
            disp(EndText);
            dispLength='%5.7f';
            disp(['InitCond=[',num2str(IC_Store(1,1),dispLength),', ',num2str(IC_Store(2,1),dispLength),', ',...
                            num2str(IC_Store(3,1),dispLength),', ',num2str(IC_Store(4,1),dispLength),', ',...
                            num2str(IC_Store(5,1),dispLength),'];',10]);
        end
    case 1
        varargout={EndReached};
        disp(EndText);
        dispLength='%5.7f';
        disp(['InitCond=[',num2str(IC_Store(1,1),dispLength),', ',num2str(IC_Store(2,1),dispLength),', ',...
                        num2str(IC_Store(3,1),dispLength),', ',num2str(IC_Store(4,1),dispLength),', ',...
                        num2str(IC_Store(5,1),dispLength),'];',10]);
    case 2
        if SimType==1
            varargout={EndReached,IC_Store(:,1)};
        else
            varargout={EndReached,max(abs(MinSlope),abs(MaxSlope))};
        end
    case 3
        varargout={IC_Store,1,1};
    case 4
        varargout={EndReached, EndText, 1, IC_Store(:,1)};
    case 5
        varargout={EndReached, EndText, T, X, TorquesP};
    otherwise
        varargout={};
end
end

