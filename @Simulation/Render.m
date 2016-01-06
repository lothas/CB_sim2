function status = Render(sim,t,X,flag)
% Renders the simulation graphics
    switch flag
        case 'init'
            t = t(1);

            if sim.Once
                % Open new figure
                if ishandle(sim.Fig) && strcmp(get(sim.Fig,'type'),'figure')
                    figure(sim.Fig);
                    if isempty(findobj(gcf,'Type','uicontrol'))
                        % Make window larger
                        set(sim.Fig,'Position', [sim.FigWidth 100,...
                            sim.FigWidth sim.FigHeight]);
                    end
                else
                    sim.Fig = figure();
                    % Make window larger
                    set(sim.Fig,'Position', [sim.FigWidth 100,...
                        sim.FigWidth sim.FigHeight]);
                end
                set(gca,'LooseInset',get(gca,'TightInset')*2)
                cla % clear previous render

                % Initialize COM tranform for "follow" mode
                sim.tCOM = hgtransform('Parent',gca);

                btn_x0 = 0.79;
                btn_y0 = 0.91;
                btn_delta = 0.0;
                if isprop(sim.Con,'s_in')
                    % Controller has higher level speed input
                    % Add controllers to simulation display
                    btn_delta = 0.05;
                    
                    uicontrol('Style', 'pushbutton', 'String', '-',...
                    'Units','normalized','FontSize',16,'FontWeight','bold',...
                    'Position', [btn_x0 btn_y0-0.05 0.035 0.035],...
                    'Callback', @sim.MinusButtonCb);
                    sim.hParam = ...
                        uicontrol('Style', 'text', 'String', 's_in = 0',...
                        'Units','normalized','FontSize',12,...
                        'Position', [btn_x0+0.04 btn_y0-0.05 0.10 0.035]);
                    uicontrol('Style', 'pushbutton', 'String', '+',...
                    'Units','normalized','FontSize',16,'FontWeight','bold',...
                        'Position', [btn_x0+0.145 btn_y0-0.05 0.035 0.035],...
                        'Callback', @sim.PlusButtonCb);
                end
                    
                % Add a 'Stop simulation' button
                uicontrol('Style', 'pushbutton', 'String', 'Stop Simulation',...
                    'Units','normalized','FontSize',12,...
                    'Position', [btn_x0 btn_y0 0.18 0.05],...
                    'Callback', @sim.StopButtonCb);
                
                % Initialize display timer
                sim.hTime = uicontrol('Style', 'text',...
                    'String', sprintf(sim.TimeStr,t,X(sim.ConCo),...
                        sim.Env.SurfSlope(sim.Mod.xS)*180/pi,sim.Mod.curSpeed),...
                    'HorizontalAlignment','left',...
                    'FontSize',12,...
                    'Units','normalized',...
                    'Position', [btn_x0+0.01 btn_y0-btn_delta-0.14 0.16 0.12],...
                    'backgroundcolor',get(gca,'color')); 
                
                % Initialize convergence display
                sim.hConv = uicontrol('Style', 'text',...
                    'String', sprintf(sim.ConvStr,1,'-'),...
                    'HorizontalAlignment','left',...
                    'FontSize',12,...
                    'Units','normalized',...
                    'Position', [btn_x0+0.01 btn_y0-btn_delta-0.2 0.16 0.06],...
                    'backgroundcolor',get(gca,'color')); 

                sim.Once = 0;

                % Render torque plots
                if sim.Con.nPulses>0
                    sim.hTorques = zeros(sim.nOuts,1);                
                    for to = 1:sim.nOuts
                        sim.hTorques(to) = line(sim.Ttime,...
                            sim.Tbase+sim.Thold(to,:)*sim.Tscale,...
                            'parent',sim.tCOM,'Color',sim.Colors{to},...
                            'LineWidth',2);
                    end
                end
            end
    end
    
    if ishandle(sim.tCOM)==0
        sim.Once = 1;
        status = Render(sim,t,X,flag);
        return
    end
    
    if ~isempty(X)
        % Check time elapsed
        sim.TimeTic = sim.TimeTic+1;
        if sim.TimeTic<sim.RSkip
            status = sim.StopSim;
            return
        else
            sim.TimeTic = 0;
        end
        
        FlMin = sim.FlMin;
        FlMax = sim.FlMax;
        if sim.Follow
            [COMx,~]=sim.Mod.GetPos(X(sim.ModCo),'COM');
            [COMy,~]=sim.Env.Surf(COMx);
            FlMin = sim.FlMin - sim.COMx0 + COMx;
            FlMax = sim.FlMax - sim.COMx0 + COMx;
            HeightMin = sim.HeightMin - sim.COMy0 + COMy;
            HeightMax = sim.HeightMax - sim.COMy0 + COMy;
            
            TCOMx = makehgtform('translate',[COMx-sim.COMx0 COMy-sim.COMy0 0]);
            set(sim.tCOM,'Matrix',TCOMx);

            axis([FlMin FlMax HeightMin HeightMax]);
        end

        % Update model render
        sim.Mod = sim.Mod.Render(X(sim.ModCo));
        % Update environment render
        sim.Env = sim.Env.Render(FlMin,FlMax);
        % Update time display
        set(sim.hTime,'string',...
            sprintf(sim.TimeStr,t(1),X(sim.ConCo),...
                sim.Env.SurfSlope(sim.Mod.xS)*180/pi,sim.Mod.curSpeed));
        % Update convergence display
        Period = find(sim.stepsSS>0,1,'first');
        if ~isempty(Period)
%             diff = norm(sim.ICstore(:,1) - sim.ICstore(:,1+Period));
            diff = min(sim.ICdiff);
            set(sim.hConv,'string',...
                sprintf(sim.ConvStr,diff,int2str(Period)),...
                    'backgroundcolor',[0.5 1 0.5]);
        else
%             diff = norm(sim.ICstore(:,1) - sim.ICstore(:,2));
            diff = min(sim.ICdiff);
            Period = find(sim.ICdiff == diff, 1, 'first');
            set(sim.hConv,'string',...
                sprintf(sim.ConvStr,diff,['(',int2str(Period),')']),...
                    'backgroundcolor',get(gca,'color'));
        end
        % Update torque display
        if sim.Con.nPulses>0
            sim.Thold(:,1:end-1) = sim.Thold(:,2:end);
            sim.Thold(:,end) = sim.Mod.Torques;
            for to = 1:sim.nOuts
                set(sim.hTorques(to),...
                    'XData',sim.Ttime,...
                    'YData',sim.Tbase + sim.Thold(to,:)*sim.Tscale);
            end
        end
    end
    status = sim.StopSim;
    drawnow
end