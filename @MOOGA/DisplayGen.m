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
                
SlopeID = 1;

% Open analysis data if it exists
FileName = ['Gen',int2str(ID),'.mat'];
if exist(FileName,'file') == 2
    % Load input file
    In = load(FileName);
    Data = In.Data;
    if ~Data.Done
        % Perform analysis and then open GUI
        disp('Performing analysis...')
        Data = GA.Analyze(Gen,ID);
    end 
else
    % Perform analysis and then open GUI
    disp('Performing analysis...')
    Data = GA.Analyze(Gen,ID);
end

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

% Add buttons
% Initial conditions
DG.menub(1) = uicontrol('Style','pushbutton',...
    'String','IC',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.88 0.08 0.05],...
    'Callback', {@Display,'IC',Data});
% Eigenvalues
DG.menub(2) = uicontrol('Style','pushbutton',...
    'String','Eigenvalues',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.82 0.08 0.05],...
    'Callback', {@Display,'EigV',Data});
% Max ZMP
DG.menub(3) = uicontrol('Style','pushbutton',...
    'String','Max ZMP',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.76 0.08 0.05],...
    'Callback', {@Display,'MZMP',Data});
% Max Torques
DG.menub(4) = uicontrol('Style','pushbutton',...
    'String','Max Torques',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.70 0.08 0.05],...
    'Callback', {@Display,'MTorques',Data});
% Power
DG.menub(5) = uicontrol('Style','pushbutton',...
    'String','Actuator Power',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.64 0.08 0.05],...
    'Callback', {@Display,'Power',Data});

% Limit cycle (states, torques, ZMP)
uicontrol('Style','Text','String','','Units','Normalized',...
    'Position', [0.02 0.623 0.08 0.001], 'BackgroundColor',[0 0 0]);
DG.menub(6) = uicontrol('Style','pushbutton',...
    'String','By slope',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.56 0.08 0.05],...
    'Callback', {@Display,'LC',Data});
% Add checkboxes
CBH = 0.03; % Height
CBS = CBH+0.01; % Spacing

DG.menub(7) = uicontrol('Style','checkbox',...
    'String','Limit Cycle','Value',1,...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.56-CBS 0.08 CBH]);
DG.menub(8) = uicontrol('Style','checkbox',...
    'String','Angles',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.56-2*CBS 0.08 CBH]);
DG.menub(9) = uicontrol('Style','checkbox',...
    'String','Ang. velocities','Value',1,...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.56-3*CBS 0.08 CBH]);
DG.menub(10) = uicontrol('Style','checkbox',...
    'String','CPG phase',...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.56-4*CBS 0.08 CBH]);
DG.menub(11) = uicontrol('Style','checkbox',...
    'String','Torques','Value',1,...
    'Units','Normalized','FontSize',ButtonFont,...
    'Position', [0.02 0.56-5*CBS 0.08 CBH]);

DG.DispArea = axes('Units','Normalized','FontSize',AxesFont,...
    'Position',[0.14 0.06 0.84 0.9]);

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
        set(hAx,'Position',[0.14 0.06 0.84 0.9]);
        % Hide legends
        hLg = findobj(gcf,'Type','axes','Tag','legend');
        set(hLg,'visible','off')
        % Hide Slope display
        hSl = FindSlopeDisp();
        if hSl
            delete(hSl)
        end
        hold on
        
        switch Type
            case 'IC'
                Ncoords = size(Data.IC,1);
                h = zeros(Ncoords,1);
                for c = 1:Ncoords
                    h(c) = plot(Data.Slopes,Data.IC(c,:),...
                        LineStyles{c},'LineWidth',LineWidth,...
                        'Color',Colors{c});
                end
                axis([min(Data.Slopes) max(Data.Slopes) ylim])
                legend(h,Legends);
            case 'EigV'
                % Separate into zones
                zID = find(diff(Data.Period(:,1))~=0);
                zID = [zID; size(Data.Period,1)];
                sID = 1;
                for z = 1:length(zID)
                    IDs = sID:zID(z);
                    sID = zID(z)+1;
                    plot(Data.Slopes(IDs),abs(Data.EigV(:,IDs)),...
                        'LineWidth',LineWidth);
                    % Add zone letter and display period number
                    zmid = mean(Data.Slopes(IDs));
                    ZString = [char(z-1+'A'),' - ',...
                        int2str(Data.Period(IDs(1),1))];
                    text(zmid,0.95,ZString,'FontSize',AxesFont,...
                        'HorizontalAlignment','center');
                    if z<length(zID)
                        % Add a vertical line
                        Lx = (Data.Slopes(sID)+Data.Slopes(sID-1))/2;
                        line([Lx Lx],[0 1],'LineWidth',LineWidth,...
                            'Color',[0 0 0])
                    end
                end
                % plot(Data.Slopes,Data.Period(:,1));
                axis([min(Data.Slopes) max(Data.Slopes) 0 1])
            case 'MZMP'
                h = zeros(3,1);
                MZMP = [Data.MZMP, Data.MZMP(:,1)-Data.MZMP(:,2)];
                for c = 1:3
                    h(c) = plot(Data.Slopes,MZMP(:,c),...
                        LineStyles{c},'LineWidth',LineWidth,...
                        'Color',Colors{c});
                end
                axis([min(Data.Slopes) max(Data.Slopes) ylim])
                legend(h,'Front','Back','Foot length');
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
                        set(hAx,'Position',[0.14 0.06 0.44 0.9]);
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

