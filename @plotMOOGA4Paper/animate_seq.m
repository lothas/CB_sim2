function [ ] = animate_seq(obj,whichGA,geneNum, GenID, Dur, timestep, filename)
%PLOT_SEQ plots the CPG output for a given sequence

GA = obj.data{1,whichGA}.GA;

Generation = geneNum;

% Animation properties
FPS = 25;
Background = [0 0 0]; % black
% Background = [1 1 1]; % white

% Render handles
AnFig = 0; tCOM = 0; COMx0 = 0; COMy0 = 0;
AnWinW = 0.4; AnWinH = 0.8; AR = 0;
FlMin = 0; FlMax = 0; HeightMin = 0; HeightMax = 0;

Sim = deepcopy(GA.Sim);
if isnumeric(Dur)
    tend = Dur;
    Sim.EndCond = 0;
    Sim.doGoNoGo = 0;
else
    switch lower(Dur)
        case 'step'
            Sim.EndCond = [1,1];
            tend = 20;
        case {'lc','limit cycle','limitcycle','limit_cycle'}
            Sim.EndCond = 2;
            tend = 'inf';
            Sim.minDiff = 1e-5;
    end
end

Sim.Graphics = 1;

Sim = GA.Gen.Decode(Sim, GA.Seqs(GenID,:,Generation));

disp(GA.Gen.seq2str(GA.Seqs(GenID,:,Generation)));

% Simulation parameters
Sim = Sim.SetTime(0,timestep,tend);

% Set internal parameters (state dimensions, events, etc)
Sim = Sim.Init();

% Some more simulation initialization
Sim.Mod.LegShift = Sim.Mod.Clearance;
Sim.Con = Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));
Sim.Con = Sim.Con.Adaptation();

% Initialize flat terrain
Sim.Env = Terrain(0,0);

% Simulate
Sim = Sim.Run();

if Sim.EndCond == 2
    Sim.Mod.xS = 0; Sim.Mod.yS = 0;
    Sim = SimStep(Sim);
end

% Clean output
Sim.Mod.Support = Sim.Mod.Right;
Sim.Mod.xS = 0; Sim.Mod.yS = 0;
T = Sim.Out.T;
X = Sim.Out.X;
Pos = Sim.Out.SuppPos;
Doubles = find(diff(T)==0);
T(Doubles) = [];
X(Doubles,:) = [];
Pos(Doubles,:) = [];

% Render animation in realtime
NPanels = floor(FPS*T(end));
dt = T(end)/NPanels;
Render(X(1,:)); % Init Render

% Foot retraction
TSteps = [T(diff(Pos(:,1))~=0); T(end)];
Nf = 226;
Shape = 1-(1-sin(0.45:0.01:2.7).^10).^10;
StepsTaken = 0;
Tf = linspace(0,TSteps(StepsTaken+1),Nf);
Sim.Mod.Clearance = 0.05*Sim.Mod.L;

OldX = X(1,:);

if ~isempty(filename)
    for p = 1:NPanels+1
        t = dt*(p-1);

        % Find the state at time t
        ID = find(abs(T-t)==min(abs(T-t)),1,'first');

        if Sim.Mod.xS~=Pos(ID,1)
            % Step taken, switch legs
            if Sim.Mod.Support==Sim.Mod.Left
                Sim.Mod.Support=Sim.Mod.Right;
            else
                if Sim.Mod.Support==Sim.Mod.Right
                    Sim.Mod.Support=Sim.Mod.Left;
                end
            end
        end

        thisX = interp1(T,X,t,'pchip');
        if norm(OldX(1:2)-thisX(1:2))>0.15
            thisX = X(ID,:);
        end
        OldX = thisX;

        % Change foot position
        Sim.Mod.xS = Pos(ID,1); Sim.Mod.yS = Pos(ID,2);

        % Retract the leg based on t
        if t>Tf(end)
            % New step
            StepsTaken = StepsTaken+1;
            Tf = linspace(TSteps(StepsTaken),TSteps(StepsTaken+1),Nf);
        end
        thisFR = interp1(Tf,Shape,t,'pchip');
        Sim.Mod.LegShift = thisFR*Sim.Mod.Clearance;
        Render(thisX);

        % Output to GIF
        frame = getframe(AnFig);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
    %     BackgroundDiff = sum(cm-repmat(Background,size(cm,1),1),2);
    %     BackgroundID = find(BackgroundDiff == min(BackgroundDiff),1,'first')-1;
        if p == 1;
            imwrite(imind,cm,filename,'gif','DelayTime',dt,'Loopcount',inf);
    %         imwrite(imind,cm,filename,'gif','DelayTime',dt,...
    %             'Loopcount',inf,'TransparentColor',BackgroundID);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',dt,'WriteMode','append');
    %         imwrite(imind,cm,filename,'gif','DelayTime',dt,...
    %             'WriteMode','append','TransparentColor',BackgroundID);
        end
    end
    disp('FINISHED RENDERING GENOME ANIMATION');

end
    function Sim = SimStep(Sim)
        Sim.EndCond = [1,2*ceil(Sim.Period(1)/2)];
        Sim.IC = Sim.IClimCyc;
        Sim.Con.Reset(Sim.IC(Sim.ConCo));
        Sim.Mod.LegShift = Sim.Mod.Clearance;
        Sim = Sim.Init();
        Sim = Sim.Run();
    end

    function Sim = SimWalk(Sim,Dur)
        if nargin<2
            Dur = 'inf';
            Sim.EndCond = 1;
        else
            if isnumeric(Dur)
                Sim.EndCond = 0;
                Sim.doGoNoGo = 0;
            else
                if strcmp(Dur,'step')
                    Dur = 10;
                    Sim.EndCond = [1,1];
                end
            end
        end
        Sim = Sim.SetTime(0,0.03,Dur);
        Sim.IC = Sim.IClimCyc;
        Sim.Con.Reset(Sim.IC(Sim.ConCo));
        Sim.Mod.LegShift = Sim.Mod.Clearance;
        close(Sim.Fig)
        Sim = Sim.Init();
        Sim = Sim.Run();
    end

    function Render(X)
        if ishandle(AnFig) && strcmp(get(AnFig,'type'),'figure')
            % Render was already initialized, just update
            figure(AnFig);
        else
            % Close Sim figure
            if ishandle(Sim.Fig)
                close(Sim.Fig)
            end
            
            % Open new figure
            AnFig = figure();
            set(AnFig,'units','normalized',...
                'Position', [0.1 0.1 AnWinW AnWinH],...
                'Color',Background);
            set(gca,'Position',[0 0 1 1],'Visible','off')
%             patch([-100 -100 100 100],[-100 100 100 -100],...
%                 -10*[1 1 1 1],Background);
            
            % Init window size params
            scrsz = get(0, 'ScreenSize');
            if scrsz(3)>2*scrsz(4) % isunix()
                % If 2 screens are used in Linux
                scrsz(3) = scrsz(3)/2;
            end
            AR = AnWinW*scrsz(3)/AnWinH/scrsz(4);

            % Initialize COM tranform for "follow" mode
            tCOM = hgtransform('Parent',gca);

%             % Render torque plots
%             if sim.Con.nPulses>0
%                 sim.hTorques = zeros(sim.nOuts,1);                
%                 for to = 1:sim.nOuts
%                     sim.hTorques(to) = line(sim.Ttime,...
%                         sim.Tbase+sim.Thold(to,:)*sim.Tscale,...
%                         'parent',sim.tCOM,'Color',sim.Colors{to},...
%                         'LineWidth',2);
%                 end
%             end
        end

        if ~isempty(X)
            % Set world view
            [COMx,~] = Sim.Mod.GetPos(X(Sim.ModCo),'COM');
            [COMy,~] = Sim.Env.Surf(COMx);
            FlMin = COMx - 0.75*AR*Sim.Mod.L;
            FlMax = COMx + 0.75*AR*Sim.Mod.L;
            HeightMin = COMy - 0.25/AR*Sim.Mod.L;
            HeightMax = COMy + 1.25/AR*Sim.Mod.L;

            TCOMx = makehgtform('translate',...
                [COMx - COMx0, COMy - COMy0, 0]);
            set(tCOM,'Matrix',TCOMx);
            axis([FlMin FlMax HeightMin HeightMax]);
            
            % Update model render
            Sim.Mod = Sim.Mod.Render(X(Sim.ModCo));
            
            % Update environment render
            Sim.Env = Sim.Env.Render(FlMin,FlMax);
            
        end
        drawnow
    end



end