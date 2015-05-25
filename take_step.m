function [X, State,failed] = take_step(Sim,IC,Slope,simT,foot_ext)
    failed = 0;
    if length(IC) == 4
        gamm = Sim.Env.incline;
        ic_temp = [IC(1)+gamm;-IC(1)+gamm;IC(2);IC(3);IC(4)];
        IC = ic_temp;
    end

    if nargin<4
        simT = {0,0.03,'inf'};
        if nargin<5
            foot_ext = 0; % foot retraction is dictated internally
        end
    end
    
    Sim.IC = IC;
    Sim.Graphics = 1;
    if strcmp(simT{3},'inf')
        Sim.EndCond = [1,1];
    else
        Sim.EndCond = 0;
    end
    % Simulation parameters
    Sim = Sim.SetTime(simT{1},simT{2},simT{3});

    % Set internal parameters (state dimensions, events, etc)
    Sim = Sim.Init();

    % Some more simulation initialization
    if ~foot_ext
        Sim.Mod.LegShift = Sim.Mod.Clearance;
    end
    % Sim.Con = Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));
    Sim.Con.lastPhi = Slope;
    Sim.Con = Sim.Con.HandleExtFB(Sim.IC(Sim.ModCo),Sim.IC(Sim.ConCo),Slope);
    
    Sim.Con = Sim.Con.Reset(Sim.IC(Sim.ConCo));
    

    % Simulate
    Sim = Sim.Run();
    if Sim.Out.Type ~= 4
        failed = 1;
    end
    State=Sim.Out;
    x = Sim.ICstore(:,1);
%   X = [(x(1)-x(2))/2;x(3);x(4);x(5)];
    X = [x(1);x(2);x(3);x(4);x(5)];
end