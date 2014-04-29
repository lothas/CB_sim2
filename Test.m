function Test()
close all

% Simulation objects
Mod = CompassBiped();
% Con = Controller();
Env = Terrain(1,0);

% Simulation parameters
tstep = 0.01;
tspan = 0:tstep:10;

% Render parameters
Once = 1;
Follow = 1;
tCOM = 0; COMx0 = 0; COMy0 = 0;
hTime = 0;
FlMin = -2; FlMax = 2;

% Initialize simulation
Mod.LegShift = Mod.Clearance;

% Simulate
options=odeset('MaxStep',tstep/10,'RelTol',.5e-7,'AbsTol',.5e-8, 'OutputFcn', @Render, 'Events', @Events);
[TTemp,XTemp,TE,YE,IE]=ode45(@Mod.Derivative,tspan,[0.1, -0.1, -0.4, -0.25],options);

while TTemp(end)<tspan(end)
    % Deal with it
    [Mod,Xa] = Mod.HandleEvent(IE,XTemp(:,end),TTemp(end));
    
    % Set new initial conditions
    IC = Xa;
    
    % Continue simulation
    [TTemp,XTemp,TE,YE,IE]=ode45(@Mod.Derivative,tspan,IC,options);
end

    function status = Render(t,X,flag)
        switch flag
            case 'init'
                t = t(1);
                
                if Once
                    figure()
                    [COMx0,COMy0]=Mod.GetPos(X,'COM');
                    tCOM = hgtransform('Parent',gca);
                    hTime = text(0.9,2.1,['t = ',num2str(t)],...
                        'HorizontalAlignment','left','Parent',tCOM);
                    Once = 0;
                end
        end
        if ~isempty(X)
            if Follow
                [COMx,~]=Mod.GetPos(X,'COM');
                [COMy,~]=Env.Surf(COMx);
                FlMin=COMx-1.25*Mod.L;
                FlMax=COMx+1.25*Mod.L;
                HeightMin=COMy-0.5*Mod.L;
                HeightMax=COMy+1.5*Mod.L;
                
                TCOMx = makehgtform('translate',[COMx-COMx0 COMy-COMy0 0]);
                set(tCOM,'Matrix',TCOMx);
                
                axis([FlMin FlMax HeightMin HeightMax]);
            end
            
            % Update model render
            Mod = Mod.Render(X);
            % Update environment render
            Env = Env.Render(FlMin,FlMax);
            % Update time display
            set(hTime,'string',['t = ',num2str(t)]);     
        end
        status = 0;
        drawnow
        
    end

    function [value isterminal direction] = Events(t, X)
        [mv mit md] = Mod.Events(X, Env);
        value = mv;
        isterminal = mit;
        direction = md;
    end
end