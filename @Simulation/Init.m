function [ sim ] = Init( sim )
    sim.CurSpeed='Computing...';

    % Set total number of simulation events
    sim.TotEvents=sim.Mod.NumEvents+sim.Con.NumEvents;

    % Set simulation time span
    sim.tspan=sim.tstart:sim.tstep:sim.tend;

%     % Call oscillator fire event
%     Events=zeros(CPG.NumEvents);
%     Events(1)=1;
%     CPG=CPG.HandleEvents(Events, tstart);
% 
%     % Set display parameters
%     global sim.StopSim;
%     StopSim=0;

%     % Scene parameters
%     if nargin<3
%         Graphics=1;
%     end
%     Follow=1;

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

