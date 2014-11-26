function Plots4Paper()
% Plotting parameters
InputFile = 'GA_11_24_13_18.mat';
SelectedGen = 4;

AutoLC = 0; % Set to 1 to plot the LC for the steepest up/downward slope
% Otherwise set the desired values here (in degrees):
DSLC = -8.35;
USLC = 7.0;    

DoPlots = [1,... Initial Conditions
           2,... Torques plot
           3,... Limit cycles
           4,... Eigenvalues vs slope
           5]; % Eigenvalues locus
       
% Plots format
AxesFont = 16;
LabelFont = 18;
TitleFont = 20;
LineWidth = 4;
LineStyles = {'-','--',':','-.','-*','-o'};
Markers = {'+','o','d','^','v'};
Colors = {[0 0 1],[1 0 0],[0 0.7 0.8],[0 0.7 0],[0.8 0 0.8]};
Legends = {'\theta_1','\theta_2','d\theta_1/dt',...
                    'd\theta_2/dt','\phi_C_P_G'};

% Close open figures
close all;

% Load data
if ~exist('GA','var')
    In = load(InputFile);
    GA = In.GA;
end
Data = GA.Analyze(GA.Progress,SelectedGen);

if AutoLC == 1
    DSLC = 1;
    USLC = length(Data.Slopes);
else
    DSLC = find(Data.Slopes >= DSLC,1,'first');
    USLC = find(Data.Slopes <= USLC,1,'last');
end
    
% Show genome results
% GA.DisplayGen(SelectedGen);

%% %%%%%%%%%%%%% Begin plots %%%%%%%%%%%%% %%

% Number of system coordinates
Ncoords = size(Data.IC,1)/max(Data.Period(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IC plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(1,DoPlots)
    % Plot angles, angular velocities and CPG phase in 3 subplots
    figure('units','normalized','Position',[0.3,0.03,0.3,0.8]);
    spheight = [4,4,2]; % Define subplot height
    Nsp = sum(spheight);

    lh = DoSubplot(1,[1,2]);
    set(gca,'XTickLabel',[]);
    ylabel('Angles [rad]','FontSize',LabelFont);
    set(lh,'Location','SouthEast');

    lh = DoSubplot(2,[3,4]);
    set(gca,'XTickLabel',[]);
    ylabel('Ang. velocities [rad/sec]','FontSize',LabelFont);
    set(lh,'Location','NorthWest');

    DoSubplot(3,5);
    xlabel('Slope [deg]','FontSize',LabelFont);
    ylabel('CPG phase','FontSize',LabelFont);
end

    function lh = DoSubplot(sp,coords)
        spID0 = sum(spheight(1:sp-1))+1;
        spIDs = spID0:spID0+spheight(sp)-1;
        subplot(Nsp,1,spIDs);
        hold on
        h = zeros(Ncoords,1);
        for p = 1:length(Data.Zones)
            pZone = Data.Zones{p};
            for z = 1:length(pZone)
                zIDs = pZone{z};
                for c = coords
                    StCoord = (p-1)*Ncoords+c;
                    h(c) = plot(Data.Slopes(zIDs),Data.IC(StCoord,zIDs),...
                        LineStyles{c},'LineWidth',LineWidth,...
                        'Color',Colors{c});
                end
            end
        end
        axis([min(Data.Slopes) max(Data.Slopes) ylim])
        lh = legend(h(coords),Legends(coords),'FontSize',AxesFont);
        set(gca,'FontSize',AxesFont);
    end

%%%%%%%%%%%%%%%%% CPG phase / Torque signal plot %%%%%%%%%%%%%%%%%
if ismember(2,DoPlots)
    % Plot the CPG phase and torque signal for a specific slope

    SlopeID = find(Data.Slopes == 0, 1, 'first');
%     SlopeID = length(Data.Slopes);

    T = Data.LCt{SlopeID};
    CPGphi = Data.LCx{SlopeID}(:,5);
    Torques = Data.LCtorques{SlopeID};

    ShiftKind = 2;
    if ShiftKind == 1
        % Shift the signal so it doesn't necessarily start right after impact
        Shift = 0.25; % Percentage of period
        ShiftN = round(Shift*length(T));
    else
        % Shift the signal so it fits a CPG cycle
        ShiftN = find(CPGphi == 0, 1, 'first');
    end

    % Find new impact time
    DeltaT = T(end) - T(ShiftN);
    ImpN = find(T >= DeltaT, 1, 'first');

    % T = [T(ShiftN:end); T(1:ShiftN-1)];
    Torques = [Torques(ShiftN:end,:); Torques(1:ShiftN-1,:)];
    CPGphi = [CPGphi(ShiftN:end,:); CPGphi(1:ShiftN-1,:)];

    figure('units','normalized','Position',[0.25,0.15,0.33,0.6]);
%     subplot(3,1,[1 2]);
    hold on
    plot(CPGphi,Torques(:,1),...
        LineStyles{2},'LineWidth',LineWidth,'Color',Colors{1});
    plot(CPGphi,Torques(:,2),...
        LineStyles{1},'LineWidth',LineWidth,'Color',Colors{2});
    hl = legend('Ankle','Hip');
    set(hl,'FontSize',LabelFont);
%     set(gca,'XTickLabel',[]);
    set(gca,'FontSize',AxesFont);
    xlabel('CPG phase','FontSize',LabelFont);
    ylabel('Torques [Nm]','FontSize',LabelFont);
    axis([CPGphi(1) CPGphi(end) ylim])

    % Add impact time
    MaxT = max(max(Torques));
    text(CPGphi(ImpN),0.3*MaxT,'Impact','rotation',90,'FontSize',LabelFont);
    DrawArrow([CPGphi(ImpN),0.3*MaxT],[CPGphi(ImpN),0.03*MaxT]);

%     subplot(3,1,3);
%     hold on
%     plot(CPGphi,CPGphi,...
%         LineStyles{3},'LineWidth',LineWidth,'Color',Colors{3});
%     set(gca,'FontSize',AxesFont);
%     xlabel('Time [sec]','FontSize',LabelFont);
%     ylabel('CPG Phase','FontSize',LabelFont);
%     axis([CPGphi(1) CPGphi(end) ylim])
end

    function DrawArrow(Tail, Head)
        Dir=Head-Tail;
        Dir=Dir/norm(Dir);
        AR = ylim/xlim;
        TrDir=[-Dir(2) Dir(1)]/AR; % transverse direction
        hLength=0.15*norm(Head-Tail);

        % Plot arrow
        Points=[Head;
                Head-hLength*Dir-hLength/4*TrDir;
                Head-hLength*Dir+hLength/4*TrDir;
                Tail];
        line(Points([1 4],1),[Points(2,2) Points(4,2)*0.9],...
            'LineWidth',3,'Color',[0 0 0]);
        patch(Points(1:3,1), Points(1:3,2),[0,0,0],'EdgeColor',[0,0,0]);
    end

%%%%%%%%%%%%%%%%%%%%%%%%% Limit cycle plot %%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(3,DoPlots)
    % Plot the limit cycle for a level slope and steepest up/downwards slope
    LCEqual = 0;
    figure('units','normalized','Position',[0.25,0.15,0.8,0.55]);
    subplot(1,9,1:3);
    PlotLC(DSLC);
    ylabel('Ang. velocity [rad/sec]','FontSize',LabelFont);

    subplot(1,9,4:6);
    PlotLC(find(Data.Slopes == 0, 1,'first'));
    xlabel('Leg Angle [rad]','FontSize',LabelFont);

    subplot(1,9,7:9);
    PlotLC(USLC);
end

    function PlotLC(Slope)
        th1 = Data.LCx{Slope}(:,1);
        th2 = Data.LCx{Slope}(:,2);
        th1t = Data.LCx{Slope}(:,3);
        th2t = Data.LCx{Slope}(:,4);
        
        % Extend plots to connect to the next
        th1e = [th1;th2(1)];
        th2e = [th2;th1(1)];
        th1te = [th1t;th2t(1)];
        th2te = [th2t;th1t(1)];
        
        hold on
        plot(th1e,th1te,LineStyles{1},'LineWidth',LineWidth,...
            'Color',Colors{1});
        plot(th2e,th2te,LineStyles{2},'LineWidth',LineWidth,...
            'Color',Colors{2});
        set(gca,'FontSize',AxesFont);
        title(sprintf('Slope: %.2f\\circ',Data.Slopes(Slope)),...
            'FontSize',LabelFont);
        if LCEqual
            axis([-0.4 0.5 -3 3])
        end
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%% Eigenvalues plot %%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(4,DoPlots)
    % Plot the Poincare map eigenvalues over the range of slopes
    figure('units','normalized','Position',[0.2,0.1,0.5,0.65]);
    % Separate into zones
    zID = find(diff(Data.Period(DSLC:USLC,1))~=0);
    zID = [zID; USLC];
    sID = DSLC;
    ZoneLetter = 0;
    for z = 1:length(zID)
        IDs = sID:zID(z);
        sID = zID(z)+1;
        if length(IDs)<2
            continue
        end
        plot(Data.Slopes(IDs),abs(Data.EigV(:,IDs)),...
            'LineWidth',LineWidth);
        if length(zID)>1
            % Add zone letter and display period number
            zmid = mean(Data.Slopes(IDs));
            ZString = [char(ZoneLetter+'A'),' - ',...
                int2str(Data.Period(IDs(1),1))];
            ZoneLetter = ZoneLetter+1;
            text(zmid,0.95,ZString,'FontSize',AxesFont,...
                'HorizontalAlignment','center');
        end
        if z<length(zID)
            % Add a vertical line
            Lx = (Data.Slopes(sID)+Data.Slopes(sID-1))/2;
            line([Lx Lx],[0 1],'LineWidth',LineWidth,...
                'Color',[0 0 0])
        end
    end
    % plot(Data.Slopes,Data.Period(:,1));
    axis([Data.Slopes(DSLC) Data.Slopes(USLC) 0 1])
    set(gca,'FontSize',AxesFont);
    xlabel('Slope [deg]')
    ylabel('|\lambda_i|')
end


%%%%%%%%%%%%%%%%%%%%% Eigenvalues locus plot %%%%%%%%%%%%%%%%%%%%%
if ismember(5,DoPlots)
    % Plot the Poincare map eigenvalues over the range of slopes
end




end