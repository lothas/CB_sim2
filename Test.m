function Test()
close all

% Simulation objects
Mod = CompassBiped();
Con = Controller();
Env = Terrain();

% Set up the compass biped model
% Mod = Mod.Set('damp',5,'yS',2.5);

% Set up the terrain
Env = Env.Set('Type','inc','start_slope',-1);

% Set up the controller
Con = Con.Set('omega0',1.106512566,'P_LegE',0.65,'FBType',0);
Con = Con.ClearTorques();
Con = Con.AddPulse('joint',1,'amp',1,'offset',1,'dur',1);
Con = Con.AddPulse('joint',1,'amp',1,'offset',1,'dur',1);

% Simulation parameters
tstep = 0.003;
tspan = 0:tstep:10;
IC = [0.13, -0.1, -0.4, -0.25];

% Render parameters
Fig = 0; Once = 1; AR = 1;
Follow = 1;
tCOM = 0; COMx0 = 0; COMy0 = 0;
hTime = 0;
FlMin = -2; FlMax = 2;

% Initialize simulation
StopSim = 0;
Mod.LegShift = Mod.Clearance;

% Simulate
options=odeset('MaxStep',tstep/10,'RelTol',.5e-7,'AbsTol',.5e-8, 'OutputFcn', @Render, 'Events', @Events);
[TTemp,XTemp,TE,YE,IE]=ode45(@Mod.Derivative,tspan,IC,options);

while TTemp(end)<tspan(end) && StopSim == 0
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
                    % Open new figure
                    Fig = figure();
                    
                    % Make window larger
                    scrsz = get(0, 'ScreenSize');
                    if isunix()
                        % If 2 screens are used in Linux
                        scrsz(3) = scrsz(3)/2;
                    end
                    FigWidth = scrsz(3)-250;
                    FigHeight = scrsz(4)-250;
                    AR = FigWidth/FigHeight;
                    set(Fig,'Position', [100 100 FigWidth FigHeight]);
                    set(gca,'LooseInset',get(gca,'TightInset')*2)
                    
                    % Initialize COM tranform for "follow" mode
                    [COMx0,COMy0]=Mod.GetPos(X,'COM');
                    tCOM = hgtransform('Parent',gca);
                    
                    % Initialize display timer
                    hTime = uicontrol('Style', 'text',...
                        'String', ['t = ',num2str(t)],...
                        'HorizontalAlignment','left',...
                        'FontSize',16,...
                        'Units','normalized',...
                        'Position', [0.88 0.85 0.1 0.03],...
                        'backgroundcolor',get(gca,'color')); 
                    
                    % Add a 'Stop simulation' button
                    uicontrol('Style', 'pushbutton', 'String', 'Stop Simulation',...
                        'Units','normalized','FontSize',16,...
                        'Position', [0.87 0.9 0.1 0.05],...
                        'Callback', @StopButtonCb);
                    Once = 0;
                end
        end
        if ~isempty(X)
            if Follow
                [COMx,~]=Mod.GetPos(X,'COM');
                [COMy,~]=Env.Surf(COMx);
                FlMin=COMx-1.25*AR*Mod.L;
                FlMax=COMx+1.25*AR*Mod.L;
                HeightMin=COMy-1/AR*Mod.L;
                HeightMax=COMy+4/AR*Mod.L;
                
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
        status = StopSim;
        drawnow
        
    end

    function [value isterminal direction] = Events(t, X)
        [mv mit md] = Mod.Events(X, Env);
        value = mv;
        isterminal = mit;
        direction = md;
    end

    function StopButtonCb(hObject, eventdata, handles) %#ok<INUSD>
        if StopSim == 0
            StopSim=1;
            set(hObject,'String','Close Window');
        else
            close(Fig)
        end
    end  % StopButtonCallback
end