function [ ] = DisplayGen(GA,ID,Gen)
%DISPLAYGEN Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
    Gen = GA.Progress;
end

TitleFont = 16;
ButtonFont = 12;
AxesFont = 12;

LineWidth = 2;
LineStyles = {'-','--',':','-.','-*','-o'};
Markers = {'+','o','d','^','v'};
Colors = {[0 0 1],[1 0 0],[0 0.7 0.8],[0 0.7 0],[0.8 0 0.8]};
Legends = {'\theta_1','\theta_2','d\theta_1/dt',...
                    'd\theta_2/dt','\phi_C_P_G'};
                

DefAxesArea = [0.14 0.06 0.84 0.9];
SlopeID = 1;

% Open analysis data if it exists
Data = GA.Analyze(Gen,ID);

DG.fig = figure();
set(DG.fig,'Toolbar','figure')

% Init window size params
scrsz = get(0, 'ScreenSize');
if scrsz(3)>2*scrsz(4) % isunix()
    % If 2 screens are used in Linux
    scrsz(3) = scrsz(3)/2;
end

% Make window larger
set(DG.fig,'Position', [100, 100, scrsz(3)-250, scrsz(4)-250]);

% Add menu box
uicontrol('Style','Text','String','Display:',...
    'Units','Normalized','FontSize',TitleFont,...
    'Position', [0.01 0.02 0.1 0.96]);

Btn0 = 0.89;
BtnH = 0.05; % Button height
BtnS = BtnH + 0.01; % Button spacing
Btn1 = Btn0 - 6*BtnS - 0.01;

% Add buttons
% Initial conditions
DG.menub(1) = uicontrol('Style','pushbutton',...
    'String','IC',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn0 0.08 BtnH],...
    'Callback', {@Display,'IC',Data});
% Eigenvalues
DG.menub(2) = uicontrol('Style','pushbutton',...
    'String','Eigenvalues',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn0-BtnS 0.08 BtnH],...
    'Callback', {@Display,'EigV',Data});
% Max ZMP
DG.menub(3) = uicontrol('Style','pushbutton',...
    'String','Max ZMP',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn0-2*BtnS 0.08 BtnH],...
    'Callback', {@Display,'MZMP',Data});
% Foot design
DG.menub(4) = uicontrol('Style','pushbutton',...
    'String','Foot design',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn0-3*BtnS 0.08 BtnH],...
    'Callback', {@Display,'FD',Data});
% Max Torques
DG.menub(5) = uicontrol('Style','pushbutton',...
    'String','Max Torques',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn0-4*BtnS 0.08 BtnH],...
    'Callback', {@Display,'MTorques',Data});
% Power
DG.menub(6) = uicontrol('Style','pushbutton',...
    'String','Actuator Power',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn0-5*BtnS 0.08 BtnH],...
    'Callback', {@Display,'Power',Data});

% Limit cycle (states, torques, ZMP)
uicontrol('Style','Text','String','','Units','Normalized',...
    'Position', [0.02 Btn1+0.0575 0.08 0.001], 'BackgroundColor',[0 0 0]);
DG.menub(7) = uicontrol('Style','pushbutton',...
    'String','By slope',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn1 0.08 BtnH],...
    'Callback', {@Display,'LC',Data});

% Add checkboxes
CBH = 0.03; % Height
CBS = CBH+0.01; % Spacing
Btn2 = Btn1 - CBS;

DG.menub(8) = uicontrol('Style','checkbox',...
    'String','Limit Cycle','Value',1,...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn2 0.08 CBH]);
DG.menub(9) = uicontrol('Style','checkbox',...
    'String','Angles',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn2-CBS 0.08 CBH]);
DG.menub(10) = uicontrol('Style','checkbox',...
    'String','Ang. velocities','Value',1,...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn2-2*CBS 0.08 CBH]);
DG.menub(11) = uicontrol('Style','checkbox',...
    'String','CPG phase',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn2-3*CBS 0.08 CBH]);
DG.menub(12) = uicontrol('Style','checkbox',...
    'String','Torques','Value',1,...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn2-4*CBS 0.08 CBH]);
DG.menub(13) = uicontrol('Style','checkbox',...
    'String','ZMP','Value',0,...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 Btn2-5*CBS 0.08 CBH]);

DG.DispArea = axes('Units','Normalized','FontSize',AxesFont,...
    'Position',DefAxesArea);

    function Display(hObject, eventdata, Type, Data) %#ok<INUSL>
        % Get the axes handle(s)
        hAx = findobj(gcf,'Type','axes','Tag','');
        
        % Clean the display before drawing
        if length(hAx)>1
            delete(hAx(2:end));
            hAx = hAx(1);
        end
        cla
        axis auto
        set(hAx,'Position',DefAxesArea);
        % Hide legends
        hLg = findobj(gcf,'Type','axes','Tag','legend');
        set(hLg,'visible','off')
        % Hide Slope display
        hSl = FindSlopeDisp();
        if ishandle(hSl) && hSl~=0
            delete(hSl)
        end
        hold on
        
        switch Type
            case 'IC'
                PlotIC(Data.Slopes,Data.IC);
            case 'EigV'
                PlotEigV(Data.Slopes,Data.EigV);
            case 'MZMP'
                PlotMaxZMP(Data.Slopes,Data.MZMP);
%                 h = zeros(3,1);
%                 MZMP = [Data.MZMP, Data.MZMP(:,1)-Data.MZMP(:,2)];
%                 for c = 1:3
%                     h(c) = plot(Data.Slopes,MZMP(:,c),...
%                         LineStyles{c},'LineWidth',LineWidth,...
%                         'Color',Colors{c});
%                 end
%                 axis([min(Data.Slopes) max(Data.Slopes) ylim])
%                 legend(h,'Front','Back','Foot length');
            case 'FD'
%                 % Make space for sizing buttons
%                 set(hAx,'Position',[DefAxesArea(1:3) 0.85]);
%                 
%                 SlopeStr = ['Slope: ',num2str(Data.Slopes(SlopeID),'%.2f')];
%                 uicontrol('Style','Text','String',SlopeStr,...
%                     'Units','Normalized','FontSize',TitleFont,...
%                     'Position', [xS 0.92 xL 0.04]);
                prompt = {'Enter ankle to toe distance:',...
                          'Enter ankle to heel distance:'...
                          'OR enter total foot length'};
                dlg_title = 'Foot design';
                num_lines = 1;
                def = {'0.15','0.15',''};
                answer = inputdlg(prompt,dlg_title,num_lines,def);
                Sizing = cellfun(@str2num,answer,'UniformOutput',0);
                Sizing = cellfun(@abs,Sizing,'UniformOutput',0);
                try
                    if isempty(Sizing{3})
                        PlotMaxZMP(Data.Slopes,Data.MZMP);
                        % Draw toe line
                        line([min(Data.Slopes) max(Data.Slopes)],...
                            [Sizing{1} Sizing{1}],...
                            'LineWidth',LineWidth,'Color',[0 0 0]);
                        % Draw heel line
                        line([min(Data.Slopes) max(Data.Slopes)],...
                            [-Sizing{2} -Sizing{2}],...
                            'LineWidth',LineWidth,'Color',[0 0 0]);
                        % Find intersections
                        mSlope = fzero(@FuncIntrsct,0,[],...
                            Data.Slopes,Data.MZMP(:,1)-Sizing{1});
                        MSlope = fzero(@FuncIntrsct,0,[],...
                            Data.Slopes,Data.MZMP(:,2)+Sizing{2});
                        text(mSlope,Sizing{1}+0.02,...
                            ['(',num2str(mSlope),',',...
                                 num2str(Sizing{1}),')']);
                        text(MSlope,-Sizing{2}+0.02,...
                            ['(',num2str(MSlope),',',...
                                 num2str(-Sizing{2}),')']);
                    else
                        % We'll start checking with the longest heel
                        if Sizing{3}>max(abs(Data.MZMP(:,2)))
                            % Don't start toe from 0
                            T0 = Sizing{3}-max(abs(Data.MZMP(:,2)));
                            H0 = max(abs(Data.MZMP(:,2)));
                        else
                            T0 = 0;
                            H0 = Sizing{3};
                        end
                        if Sizing{3}>max(Data.MZMP(:,1))-T0
                            aMax = (max(Data.MZMP(:,1))-T0)/Sizing{3};
                        else
                            aMax = 0.99;
                        end
                        
                        % Check possible slopes
                        Npoints = 100;
                        dL = linspace(0.01,aMax,Npoints)*(Sizing{3}-T0);
                        Toe = zeros(Npoints,1);
                        Heel = zeros(Npoints,1);
                        Slopes = zeros(Npoints,3);
                        for i = 1:Npoints
                            Toe(i) = T0+dL(i);
                            Heel(i) = H0-dL(i);
                            % Negative slope
                            Slopes(i,1) = fzero(@FuncIntrsct,0,[],...
                                Data.Slopes,Data.MZMP(:,1)-Toe(i));
                            % Positive slope
                            Slopes(i,2) = fzero(@FuncIntrsct,0,[],...
                                Data.Slopes,Data.MZMP(:,2)+Heel(i));
                            % Range
                            Slopes(i,3) = Slopes(i,2) - Slopes(i,1);
%                             disp(['Foot size: ',num2str(Toe(i)+Heel(i),'%.2f'),...
%                                 ' - Min/Max slope: ',...
%                                 num2str(Slopes(i,1),'%.2f'),' / '...
%                                 num2str(Slopes(i,2),'%.2f')]);
                        end
                        
                        % Display results
                        % Make space for axes labels
                        Area = DefAxesArea + [0 0.05 0 -0.05];
                        
                        delete(hAx)
                        hAx = axes('Position',Area,...
                                      'XAxisLocation','bottom');                        
                        plot(Toe,abs(Slopes));
                        axis([min(Toe) max(Toe) ylim])
                        xlabel('Ankle to Toe [m]')
                        legend('Ankle to Toe slopes',...
                            'Ankle to Heel slopes',...
                            'Range');
                        disp(['Foot size(s): ',num2str(unique(Toe+Heel))])
                    end
                catch err
                    disp('ERROR: Wrong input');
                    disp(err)
                end                    
                
            case 'MTorques'
                h = zeros(2,1);
                for c = 1:2
                    h(c) = plot(Data.Slopes,Data.MTorques(:,c),...
                        LineStyles{c},'LineWidth',LineWidth,...
                        'Color',Colors{c});
                end
                axis([min(Data.Slopes) max(Data.Slopes) ylim])
                legend(h,'Ankle','Hip');
            case 'Power'
                h = zeros(2,1);
                for c = 1:2
                    h(c) = plot(Data.Slopes,Data.Power(:,c),...
                        LineStyles{c},'LineWidth',LineWidth,...
                        'Color',Colors{c});
                end
                axis([min(Data.Slopes) max(Data.Slopes) ylim])
                legend(h,'Ankle','Hip');
            case 'LC'
                % Get selected checkboxes
                checks = GetSelected();
                
                xS = 0.14; % axes x start
                xL = 0.84; % axes length
                
                LCID = find(ismember(checks,'Limit Cycle'),1);
                if ~isempty(LCID)
                    if length(checks)>1
                        % Plot limit cycle in smaller area
                        set(hAx,'Position',[DefAxesArea(1:2) 0.44 0.9]);
                        xS = 0.62;
                        xL = 0.34;
                        checks(LCID) = [];
                    end
                    PlotLC(Data.LCx{SlopeID});
                else
                    delete(gca)
                end
                
                NPlots = length(checks);
                AxGap = 0.04;
                AxHeight = (0.9-NPlots*AxGap)/NPlots;
                AxGH = AxGap+AxHeight;
                
                for p = 1:NPlots
                    axes('Units','Normalized','FontSize',AxesFont,...
                        'Position',[xS 0.06+(p-1)*AxGH xL AxHeight]);
                    if strcmp(checks{p},'Angles')
                        % Plot angles
                        PlotStates(Data.LCt{SlopeID},Data.LCx{SlopeID},[1 2]);
                    end

                    if strcmp(checks{p},'Ang. velocities')
                        % Plot angular velocities
                        PlotStates(Data.LCt{SlopeID},Data.LCx{SlopeID},[3 4]);
                    end
                    
                    if strcmp(checks{p},'CPG phase')
                        % Plot CPG phase
                        PlotStates(Data.LCt{SlopeID},Data.LCx{SlopeID},5);
                    end

                    if strcmp(checks{p},'Torques')
                        % Plot Torques
                        PlotTorques(Data.LCt{SlopeID},Data.LCtorques{SlopeID});
                    end
                    
                    if strcmp(checks{p},'ZMP')
                        % Plot Torques
                        PlotZMP(Data.LCt{SlopeID},Data.LCZMP{SlopeID});
                    end
                    
                    if p<NPlots
                        set(gca,'XTickLabel','')
                    end
                end
                
                % Show current slope
                SlopeStr = ['Slope: ',num2str(Data.Slopes(SlopeID),'%.2f')];
                uicontrol('Style','Text','String',SlopeStr,...
                    'Units','Normalized','FontSize',TitleFont,...
                    'Position', [xS 0.92 xL 0.04]);
                
                % Set keyboard callback
                set(gcf,'KeyPressFcn',{@KeybHandle,Data})
        end
    end

    function PlotIC(Slopes,IC)
        Ncoords = size(IC,1)/max(Data.Period(:,1));
        h = zeros(Ncoords,1);
        for p = 1:length(Data.Zones)
            pZone = Data.Zones{p};
            for z = 1:length(pZone)
                Coords = pZone{z};
                for c = 1:Ncoords
                    StCoord = (p-1)*Ncoords+c;
                    h(c) = plot(Slopes(Coords),IC(StCoord,Coords),...
                        LineStyles{c},'LineWidth',LineWidth,...
                        'Color',Colors{c});
                end
            end
        end
        axis([min(Slopes) max(Slopes) ylim])
        legend(h,Legends);
    end

    function PlotEigV(Slopes,EigV)
        % Separate into zones
        zID = find(diff(Data.Period(:,1))~=0);
        zID = [zID; size(Data.Period,1)];
        sID = 1;
        ZoneLetter = 0;
        for z = 1:length(zID)
            IDs = sID:zID(z);
            sID = zID(z)+1;
            if length(IDs)<2
                continue
            end
            plot(Slopes(IDs),abs(EigV(:,IDs)),...
                'LineWidth',LineWidth);
            % Add zone letter and display period number
            zmid = mean(Slopes(IDs));
            ZString = [char(ZoneLetter+'A'),' - ',...
                int2str(Data.Period(IDs(1),1))];
            ZoneLetter = ZoneLetter+1;
            text(zmid,0.95,ZString,'FontSize',AxesFont,...
                'HorizontalAlignment','center');
            if z<length(zID)
                % Add a vertical line
                Lx = (Slopes(sID)+Slopes(sID-1))/2;
                line([Lx Lx],[0 1],'LineWidth',LineWidth,...
                    'Color',[0 0 0])
            end
        end
        % plot(Data.Slopes,Data.Period(:,1));
        axis([min(Slopes) max(Slopes) 0 1])
    end

    function PlotMaxZMP(Slopes,MZMP)
        h = zeros(3,1);
        MZMP = [MZMP, MZMP(:,1)-MZMP(:,2)];
        pZone = Data.Zones{1};
        for z = 1:length(pZone)
            Coords = pZone{z};
            for c = 1:3
                h(c) = plot(Slopes(Coords),MZMP(Coords,c),...
                    LineStyles{c},'LineWidth',LineWidth,...
                    'Color',Colors{c});
            end
        end
        axis([min(Slopes) max(Slopes) ylim])
        legend(h,'Front','Back','Foot length');
    end

    function PlotLC(X)
        X = [X; X(1,[2 1 4 3 5])];
        hold on
        plot(X(:,1),X(:,3),...
            LineStyles{1},'LineWidth',LineWidth,'Color',Colors{1});
        plot(X(:,2),X(:,4),...
            LineStyles{2},'LineWidth',LineWidth,'Color',Colors{2});
        legend('Stance','Swing');
    end

    function PlotStates(T,X,which)
        X = [X; X(1,[2 1 4 3 5])];
        T = [T; T(end)];
        hold on
        h = zeros(size(which));
        for s = 1:length(which)
            h(s) = plot(T,X(:,which(s)),...
                LineStyles{which(s)},'LineWidth',LineWidth,...
                'Color',Colors{which(s)});
        end
        legend(h,Legends(which));
    end

    function PlotTorques(T,Torques)
        hold on
        plot(T,Torques(:,1),...
            LineStyles{1},'LineWidth',LineWidth,'Color',Colors{1});
        plot(T,Torques(:,2),...
            LineStyles{2},'LineWidth',LineWidth,'Color',Colors{2});
        legend('Ankle','Hip');
    end

    function PlotZMP(T,ZMP)
        hold on
        plot(T,ZMP,...
            LineStyles{1},'LineWidth',LineWidth,'Color',Colors{1});
    end

    function h = FindSlopeDisp()
        hs = findobj(gcf,'Type','uicontrol','Style','Text');
        for i = 1:length(hs)
            if ~isempty(strfind(get(hs(i),'String'),'Slope'))
                h = hs(i);
                return
            end
        end
        % SlopeDisp object not found
        h = 0;
    end

    function val = FuncIntrsct(x,X,Y)
        if x<min(X)
            x = min(X);
        end
        if x>max(X)
            x = max(X);
        end
        val = interp1(X,Y,x);
    end

    function checks = GetSelected()
        hs = findobj(gcf,'Type','uicontrol','Style','checkbox');
        checks = {};
        for i = 1:length(hs)
            if get(hs(i),'Value')
                checks{end+1} = get(hs(i),'String'); %#ok<AGROW>
            end
        end
    end

    function KeybHandle(h_obj,evt,Data) %#ok<INUSL>
        Sections = length(Data.Slopes);
        BigJump = floor(Sections/10);
        
        if strcmp(evt.Key,'uparrow')
            SlopeID=SlopeID+1;
            if SlopeID>Sections
                SlopeID=Sections;
            end
        end
        
        if strcmp(evt.Key,'downarrow')
            SlopeID=SlopeID-1;
            if SlopeID<1
                SlopeID=1;
            end     
        end
        
        if strcmp(evt.Key,'rightarrow')
            SlopeID=SlopeID+BigJump;
            if SlopeID>Sections
                SlopeID=Sections;
            end
        end
        
        if strcmp(evt.Key,'leftarrow')
            SlopeID=SlopeID-BigJump;
            if SlopeID<1
                SlopeID=1;
            end     
        end

        Display(1, 1, 'LC', Data);
%         X=LimitCycle_X{IPPIndex};
%         Slope=Slopes(IPPIndex);
%         
%         set(IPPh1,'XData',X(:,3));
%         set(IPPh1,'YData',X(:,5));
%         set(IPPh2,'XData',X(:,4));
%         set(IPPh2,'YData',X(:,6));
%         set(IPPt1,'String',sprintf('Interactive Phase Plane - Slope=%.2f\\circ',Slope));
%         
%         if Plots(5)==1 % If interactive PM locus is active, apply change to it as well
%             if ~isfield(evt,'Called') % check if this handle was called from the other handle
%                 evt.Called=1;
%                 PMKeybHandle([],evt);
%             end
%         end
    end
end

