function Test()
close all

% Simulation objects
Mod = CompassBiped();
Con = Controller();
Env = Terrain();

% Set up the compass biped model
% Mod = Mod.Set('damp',5,'yS',2.5);

% Set up the terrain
Env = Env.Set('Type','inc','start_slope',0);

% Set up the controller
% Con = Con.ClearTorques();
% Con = Con.Set('omega0',1.2666,'P_LegE',0.5973,'FBType',0);
% Con = Con.AddPulse('joint',1,'amp',-7.3842,'offset',0.1268,'dur',0.07227);
% Con = Con.AddPulse('joint',2,'amp',5.1913,'offset',0.1665,'dur',0.0537);

% Set states
stDim = Mod.stDim + Con.stDim;
ModCo = 1:Mod.stDim; % Model coord. indices
ConCo = Mod.stDim+1:stDim; % Contr. coord. indices

% Set events
nEvents = Mod.nEvents + Con.nEvents;
ModEv = 1:Mod.nEvents; % Model events indices
ConEv = Mod.nEvents+1:nEvents; % Contr. events indices

% Simulation parameters
tstep = 0.05;
tend = 10;
tspan = 0:tstep:tend;
% IC = [0.13, -0.1, -0.4, -0.25, 0];
IC = [0., 0., 0., 0., 0.];
IC = [0.1393442, -0.1393442, -0.5933174, -0.4680616, 0.8759402];

% Render parameters
Fig = 0; Once = 1; AR = 1;
Follow = 1;
tCOM = 0; COMx0 = 0; COMy0 = 0;
hTime = 0;
FlMin = -2; FlMax = 2;
TimeStr = 't = %.2f \nOsc. = %.2f';
% Torque render
nOuts = size(Con.OutM,1);
nTsteps = 100;
Thold = zeros(nOuts,nTsteps);
Tbase = 0; Tscale = 1;
hTorques = zeros(nOuts,1);
Colors = {[1 0 0],[0 0 1],[0 1 0],[0 0 0]};

% Initialize simulation
StopSim = 0;
Mod.LegShift = Mod.Clearance;
Con.HandleEvent(1, IC(ConCo));

% Simulate
options=odeset('MaxStep',tstep/10,'RelTol',.5e-7,'AbsTol',.5e-8, 'OutputFcn', @Render, 'Events', @Events);
[TTemp,XTemp,TE,YE,IE]=ode45(@Derivative,tspan,IC,options);

while TTemp(end)<tspan(end) && StopSim == 0
    % Deal with it
    Xa = XTemp(end,:);
    for ev = 1:length(IE)
        % Is it a model event?
        ModEvID = find(IE(ev) == ModEv,1,'first');
        if ~isempty(ModEvID)
            [Mod,Xa(ModCo)] = Mod.HandleEvent(ModEvID, ...
                                    XTemp(end,ModCo),TTemp(end));
            
            % Handle event interactions
            if ModEvID == 1
                % Ground contact
                [Con, Xa] = Con.HandleExtFB(Xa);
            end
        end
        
        % Is it a controller event?
        ConEvID = find(IE(ev) == ConEv,1,'first');
        if ~isempty(ConEvID)
            [Con,Xa(ConCo)] = Con.HandleEvent(ConEvID, XTemp(end,ConCo));
            
            % Handle event interactions
            switch ConEvID
                case 1 % Neuron fired
                    Mod.LegShift = Mod.Clearance;
                case 2 % Leg extension
                    Mod.LegShift = 0;
            end 
        end
    end
    
    % Set new initial conditions
    IC = Xa;
    
    % Continue simulation
    tspan = TTemp(end):tstep:tend;
    [TTemp,XTemp,TE,YE,IE]=ode45(@Derivative,tspan,IC,options);
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
                        'String', sprintf(TimeStr,t,X(ConCo)),...
                        'HorizontalAlignment','left',...
                        'FontSize',12,...
                        'Units','normalized',...
                        'Position', [0.88 0.8 0.1 0.08],...
                        'backgroundcolor',get(gca,'color')); 
                    
                    % Add a 'Stop simulation' button
                    uicontrol('Style', 'pushbutton', 'String', 'Stop Simulation',...
                        'Units','normalized','FontSize',12,...
                        'Position', [0.87 0.9 0.1 0.05],...
                        'Callback', @StopButtonCb);
                    Once = 0;
                    
                    % Render torque plots
                    FlMin=COMx0-1.25*AR*Mod.L;
                    FlMax=COMx0+1.25*AR*Mod.L;
                    HeightMin=COMy0-1/AR*Mod.L;
                    HeightMax=COMy0+4/AR*Mod.L;
                    Ttime = linspace(FlMax*0.8,FlMax*0.95,nTsteps);
                    Tbase = (HeightMax+HeightMin)/2;
                    Tscale = 0.1*(HeightMax-HeightMin)/max(abs(Con.Amp0));
                    for to = 1:nOuts
                        hTorques(to) = line(Ttime,...
                            Tbase+Thold(to,:)*Tscale,...
                            'parent',tCOM,'Color',Colors{to},...
                            'LineWidth',2);
                    end
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
            set(hTime,'string',sprintf(TimeStr,t,X(ConCo)));
            % Update torque display
            Thold(:,1:end-1) = Thold(:,2:end);
%             Thold(:,end) = Con.NeurOutput();
            Thold(:,end) = Mod.GetTorques();
            for to = 1:nOuts
                set(hTorques(to),'YData',Tbase + Thold(to,:)*Tscale);
            end
        end
        status = StopSim;
        drawnow
        
    end

    function [Xt] = Derivative(t,X)
        Mod.Torques=Con.NeurOutput();
        
        Xt = [Mod.Derivative(t,X(ModCo));
        	  Con.Derivative(t,X(ConCo))];
    end

    function [value, isterminal, direction] = Events(t, X) %#ok<INUSL>
        value = zeros(nEvents,1);
        isterminal = ones(nEvents,1);
        direction = zeros(nEvents,1);
        
        % Call model event function
        [value(ModEv), isterminal(ModEv), direction(ModEv)] = ...
            Mod.Events(X(ModCo), Env);
        % Call controller event function
        [value(ConEv), isterminal(ConEv), direction(ConEv)] = ...
            Con.Events(X(ConCo));
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