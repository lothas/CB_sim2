function [ res ] = Run( )
    % Model Params
    MP = 0;
    
    % Controller Params
    CP = 0;
    
    % Terrain Params
    TP = 0;
    
    [Robot, CPG, Terrain] = SimInit(MP, CP, TP);
    Sim = Simulation(Robot, CPG, Terrain);

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
end

