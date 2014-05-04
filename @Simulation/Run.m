function [ X, T ] = Run( sim )
% Run the simulation until an event occurs
% Handle the event and keep running
    X = [];
    T = [];
    
    if sim.Graphics == 1
        options=odeset('MaxStep',sim.tstep/10,'RelTol',.5e-7,'AbsTol',.5e-8,...
            'OutputFcn', @sim.Render, 'Events', @sim.Events);
    else
        options=odeset('MaxStep',sim.tstep/10,'RelTol',.5e-7,'AbsTol',.5e-8,...
            'Events', @sim.Events);
    end
    
    tspan = sim.tstart:sim.tstep:sim.tend;
    [TTemp,XTemp,TE,YE,IE] = ...
        ode45(@sim.Derivative,tspan,sim.IC,options); %#ok<ASGLU>
    
    if sim.infTime == 1
        Cond1 = true;
        sim.tend = sim.tend + TTemp(end)-tspan(1);
    else
        Cond1 = TTemp(end)<sim.tend;
    end
    
    % Save state and time
    X = [X; XTemp];
    T = [T; TTemp];

    while Cond1 && sim.StopSim == 0
        % Deal with it
        Xa = XTemp(end,:);
        for ev = 1:length(IE)
            % Is it a model event?
            ModEvID = find(IE(ev) == sim.ModEv,1,'first');
            if ~isempty(ModEvID)
                [sim.Mod,Xa(sim.ModCo)] = ...
                    sim.Mod.HandleEvent(ModEvID, ...
                        XTemp(end,sim.ModCo),TTemp(end));

                % Handle event interactions
                if ModEvID == 1
                    % Ground contact
                    [sim.Con, Xa] = sim.Con.HandleExtFB(Xa);
                    
                    % Count step
                    sim.StepsTaken = sim.StepsTaken+1;
                    
                    if ~ischar(sim.Mod.curSpeed)
                        sim.TimeStr = ['t = %.2f s\nOsc.=%.3f\n',...
                         'Slope = %.2f ',char(176)','\nSpeed = %.3f m/s'];
                    end
                end
            end

            % Is it a controller event?
            ConEvID = find(IE(ev) == sim.ConEv,1,'first');
            if ~isempty(ConEvID)
                [sim.Con,Xa(sim.ConCo)] = ...
                    sim.Con.HandleEvent(ConEvID, XTemp(end,sim.ConCo));

                % Handle event interactions
                switch ConEvID
                    case 1 % Neuron fired
                        sim.Mod.LegShift = sim.Mod.Clearance;
                    case 2 % Leg extension
                        sim.Mod.LegShift = 0;
                end 
            end
        end

        % Set new initial conditions
        sim.IC = Xa;

        % Continue simulation
        tspan = TTemp(end):sim.tstep:sim.tend;
        [TTemp,XTemp,TE,YE,IE] = ...
            ode45(@sim.Derivative,tspan,sim.IC,options); %#ok<ASGLU>
        
        if sim.infTime == 1
            Cond1 = true;
            sim.tend = sim.tend + TTemp(end)-tspan(1);
        else
            Cond1 = TTemp(end)<sim.tend;
        end
        
        % Save state and time
        X = [X; XTemp]; %#ok<AGROW>
        T = [T; TTemp]; %#ok<AGROW>
    end

end

