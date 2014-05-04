function status = Render(sim,t,X,flag)
% Renders the simulation graphics
    switch flag
        case 'init'
            t = t(1);

            if sim.Once
                % Open new figure
                sim.Fig = figure();

                % Make window larger
                set(sim.Fig,'Position', [100 100 sim.FigWidth sim.FigHeight]);
                set(gca,'LooseInset',get(gca,'TightInset')*2)

                % Initialize COM tranform for "follow" mode
                sim.tCOM = hgtransform('Parent',gca);

                % Initialize display timer
                sim.hTime = uicontrol('Style', 'text',...
                    'String', sprintf(sim.TimeStr,t,X(sim.ConCo),...
                        sim.Env.SurfSlope(sim.Mod.xS)*180/pi,sim.Mod.curSpeed),...
                    'HorizontalAlignment','left',...
                    'FontSize',12,...
                    'Units','normalized',...
                    'Position', [0.88 0.76 0.1 0.12],...
                    'backgroundcolor',get(gca,'color')); 

                % Add a 'Stop simulation' button
                uicontrol('Style', 'pushbutton', 'String', 'Stop Simulation',...
                    'Units','normalized','FontSize',12,...
                    'Position', [0.87 0.9 0.1 0.05],...
                    'Callback', @sim.StopButtonCb);
                sim.Once = 0;

                % Render torque plots
                for to = 1:sim.nOuts
                    sim.hTorques(to) = line(sim.Ttime,...
                        sim.Tbase+sim.Thold(to,:)*sim.Tscale,...
                        'parent',sim.tCOM,'Color',sim.Colors{to},...
                        'LineWidth',2);
                end
            end
    end
    
    if ~isempty(X)
        if sim.Follow
            [COMx,~]=sim.Mod.GetPos(X(sim.ModCo),'COM');
            [COMy,~]=sim.Env.Surf(COMx);
            sim.FlMin=COMx-1.25*sim.AR*sim.Mod.L;
            sim.FlMax=COMx+1.25*sim.AR*sim.Mod.L;
            HeightMin=COMy-1/sim.AR*sim.Mod.L;
            HeightMax=COMy+4/sim.AR*sim.Mod.L;

            TCOMx = makehgtform('translate',[COMx-sim.COMx0 COMy-sim.COMy0 0]);
            set(sim.tCOM,'Matrix',TCOMx);

            axis([sim.FlMin sim.FlMax HeightMin HeightMax]);
        end

        % Update model render
        sim.Mod = sim.Mod.Render(X(sim.ModCo));
        % Update environment render
        sim.Env = sim.Env.Render(sim.FlMin,sim.FlMax);
        % Update time display
        set(sim.hTime,'string',...
            sprintf(sim.TimeStr,t,X(sim.ConCo),...
                sim.Env.SurfSlope(sim.Mod.xS)*180/pi,sim.Mod.curSpeed));
        % Update torque display
        if sim.nOuts>0
            sim.Thold(:,1:end-1) = sim.Thold(:,2:end);
    %             sim.Thold(:,end) = sim.Con.NeurOutput();
            sim.Thold(:,end) = sim.Mod.Torques;
            for to = 1:sim.nOuts
                set(sim.hTorques(to),...
                    'YData',sim.Tbase + sim.Thold(to,:)*sim.Tscale);
            end
        end
    end
    status = sim.StopSim;
    drawnow
end