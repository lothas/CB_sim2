function Plots()
% Version 0.4 - 14/06/2012
% close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% To get great plots for printing use these 2 lines for each figure:     %
% set(gcf,'PaperType','A2')                                              %
% print -dtiffn -r600                                                    %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Set format parameters
FontSize=18;
LineWidth=2;

% Select to display slope in radians or degrees
RadORDeg=1; % 1 for Deg, 2 for Rad
RMult=pi/180;

% Select to display leg angles from perpendicular
% to the horizon or to the surface
% 1 for symmetrical, 0 for real values
Symmetric=0;

Plots=zeros(1,19);
% Plots to display:
% Plots(1)=1;     % Initial conditions: q_1 q_2 \phi
% Plots(2)=1;     % Initial conditions: q_1dot q_2dot
% Plots(3)=1;     % Interactive Phase Plane
% Plots(4)=1;     % Poincare eigenvalues
% Plots(5)=1;     % Interactive Poincare Map
% Plots(6)=1;     % Poincare Map eigenvalues locus
% Plots(7)=1;     % Limit Cycle plots for paper
% Plots(8)=1;     % Compass Gait Model Figure
% Plots(9)=1;     % Interactive GRF plot
% Plots(10)=1;    % Energy plots
Plots(11)=1;    % Mega interactive plot
Plots(12)=1;    % ZMP and friction coeff. plot
% %%%% !!!! Following plots not yet updated !!!! %%%% %
% Plots(13)=1;     % Muscle Torques
% Plots(14)=1;     % Potential, Kinetic and Total Energy
% Plots(15)=1;     % Energy Loss by Leg Apperture
% Plots(16)=1;    % Off COM Distance
% Plots(17)=1;    % Off COM Distance by Slope and Leg Apperture
% Plots(18)=1;    % Plots for paper
% Plots(19)=1;    % Stick-figure images

% If both initial condition plots are selected define
% whether to plot them as separate figures or as subplots
InitCondSubplots=1;

% If both interactive plots are selected define
% whether to plot them as separate figures or as subplots
InterSubplots=1;

% Define whether to print PM locus as subplots
LocusSubplots=1;
% Define whether to add direction to the locus (increasing with slope)
LocusDirection=1;
% Define whether to add an index number to the direction arrows
LocusDirNumber=0;

S=load('NumericalData.mat');   % Load precomputed numerical data
Slopes=S.Slopes;
LimitCycle_Period=S.LimitCycle_Period;
LimitCycle_IC=S.LimitCycle_IC;
LimitCycle_X=S.LimitCycle_X;
LimitCycle_T=S.LimitCycle_T;
LimitCycle_Torques=S.LimitCycle_Torques;
LimitCycle_GRF=S.LimitCycle_GRF;
LimitCycle_ZMP=S.LimitCycle_ZMP;
LimitCycle_eig=S.LimitCycle_eig;
MaxZMP = S.MaxZMP;
MaxFriction = S.MaxFriction;

StepLength=S.StepLength;
dKineticEIm=S.dKineticEIm;
dKineticEFr=S.dKineticEFr;
dPotentialE=S.dPotentialE;
ControllerE=S.ControllerE;
% AbsControllerE=S.AbsControllerE;
DPSlopes=S.DPSlopes;
DPStepLength=S.DPStepLength;
DPdKineticEIm=S.DPdKineticEIm;
DPdKineticEFr=S.DPdKineticEFr;
DPdPotentialE=S.DPdPotentialE;
DPControllerE=S.DPControllerE;
% DPAbsControllerE=S.DPAbsControllerE;

GRF=S.GRF;
MeasuredGRF=S.MeasuredGRF;
MeanMeasuredGRF=S.MeanMeasuredGRF;
MaxGRFRatio=S.MaxGRFRatio;

Sections=length(Slopes);    % Obtain number of sections

IPPIndex=ceil(Sections/2);     % Set interactive phase plane plot index to middle of range
IPMIndex=ceil(Sections/2);     % Set interactive PM locus plot index to middle of range
IGRFIndex=ceil(Sections/2);     % Set interactive GRF plot index to middle of range
IMIPIndex=ceil(Sections/2);     % Set mega interactive plot index to middle of range

% Apply modifications to display units as deg or rad
if RadORDeg==1
    SlopeLabel='Slope [deg]';
end
if RadORDeg==2
    Slopes=Slopes*RMult;
    SlopeLabel='Slope [rad]';
    RMult=1;
end

% Assign sets of data to each coordinate
nSets=floor(size(LimitCycle_IC,1)/5);
switch nSets
    case 1
        Coords=(1:5)';
    case 2
        Coords=[1:5; 6:10]';
    case 3
        Coords=[1:5; 6:10; 11:15]';
    case 4
        Coords=[1:5; 6:10; 11:15; 16:20]';
end

% We know that there's only one region where there is a period
% bifurcation so we'll split the slopes into 3 regions
Split1=find(LimitCycle_Period==2,1,'first');
Split2=find(LimitCycle_Period==2,1,'last');

Region1=1:Split1-1;
Region2=Split1:Split2;
Region3=Split2+1:length(Slopes);

Colors={[0.7 0.7 0], [1 0 1], [0 0.7 0.7], [0 0.7 0], [0 0 1], [0 0 0], [0.5 1 0.5], [0.5 0 1], [0.5 0.5 1], [1 0.5 0.25]};
LineStyles={'-', '-.', '--', '-.', '--'};
        
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Initial conditions: q_1 q_2 \phi %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Plots(1)
    figure(1);
    if InitCondSubplots==1 && Plots(2)==1
        subplot(8,1,1:3);
    end
    if nSets>1
        h1=plot(Slopes,LimitCycle_IC(Coords(1,1),:)-Symmetric*repmat(Slopes*RMult,1,1),'b','LineWidth',LineWidth);
        hold on
        h2=plot(Slopes,LimitCycle_IC(Coords(2,1),:)-Symmetric*repmat(Slopes*RMult,1,1),'--r','LineWidth',LineWidth);
        for i=2:nSets
            RepIndex=find(abs(LimitCycle_IC(Coords(1,1),:)-LimitCycle_IC(Coords(1,i),:))>1e-6);
            plot(Slopes(RepIndex),LimitCycle_IC(Coords(1,i),RepIndex)-Symmetric*repmat(Slopes(RepIndex)*RMult,1,1),'b','LineWidth',LineWidth);
            RepIndex=find(abs(LimitCycle_IC(Coords(2,1),:)-LimitCycle_IC(Coords(2,i),:))>1e-6);
            plot(Slopes(RepIndex),LimitCycle_IC(Coords(2,i),RepIndex)-Symmetric*repmat(Slopes(RepIndex)*RMult,1,1),'--r','LineWidth',LineWidth);
        end
    else
        h1=plot(Slopes,LimitCycle_IC(Coords(1,:),:)-Symmetric*repmat(Slopes*RMult,nSets,1),'b','LineWidth',LineWidth);
        hold on
        h2=plot(Slopes,LimitCycle_IC(Coords(2,:),:)-Symmetric*repmat(Slopes*RMult,nSets,1),'--r','LineWidth',LineWidth);
    end
            
    if InitCondSubplots==1 && Plots(2)==1
        t1=title('Steady state initial conditions: \theta_1,\theta_2,\theta_1_t,\theta_2_t,\phi versus slope');
        set(gca,'XTickLabel','');
    else        
        t1=title('Steady state initial conditions for \theta_1,\theta_2 versus slope');
        xh=xlabel(SlopeLabel);
        set(xh,'FontSize',FontSize);
    end
    yh=ylabel('Angle [rad]');
    legend([h1(1);h2(1)],{'\theta_1';'\theta_2'},'Orientation','Horizontal','Location','SouthEast');
    axis([min(Slopes) max(Slopes) -0.5 0.4]);
    set(gca,'FontSize',FontSize-2);
    set(yh,'FontSize',FontSize);
    set(t1,'FontSize',FontSize);
    
    if InitCondSubplots==1 && Plots(2)==1
        subplot(8,1,7:8);
    else
        figure(2);
    end
    if nSets>1
        plot(Slopes,LimitCycle_IC(Coords(5,1),:),'k','LineWidth',LineWidth);
        hold on
        for i=2:nSets
            RepIndex=find(abs(LimitCycle_IC(Coords(5,1),:)-LimitCycle_IC(Coords(5,i),:))>1e-6);
            plot(Slopes(RepIndex),LimitCycle_IC(Coords(5,i),RepIndex),'k','LineWidth',LineWidth);
        end
    else
        plot(Slopes,LimitCycle_IC(Coords(5,:),:),'k','LineWidth',LineWidth);
    end
    if InitCondSubplots==0 || Plots(2)==0
        t1=title('Steady state initial conditions for CPG phase \phi versus slope');
        set(t1,'FontSize',FontSize);
    end
    xh=xlabel(SlopeLabel);
    set(xh,'FontSize',FontSize);
    yh=ylabel('CPG Phase \phi');
    axis([min(Slopes) max(Slopes) 0.65 1]);
    set(gca,'FontSize',FontSize-2);
    set(yh,'FontSize',FontSize);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Initial conditions: q_1dot q_2dot %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Plots(2)
    if InitCondSubplots==1 && Plots(1)==1
        subplot(8,1,4:6);
    else
        figure(3);
    end
    if nSets>1
        h1=plot(Slopes,LimitCycle_IC(Coords(3,1),:),'b','LineWidth',LineWidth);
        hold on
        h2=plot(Slopes,LimitCycle_IC(Coords(4,1),:),'--r','LineWidth',LineWidth);
        for i=2:nSets
            RepIndex=find(abs(LimitCycle_IC(Coords(3,1),:)-LimitCycle_IC(Coords(3,i),:))>1e-6);
            plot(Slopes(RepIndex),LimitCycle_IC(Coords(3,i),RepIndex),'b','LineWidth',LineWidth);
            RepIndex=find(abs(LimitCycle_IC(Coords(4,1),:)-LimitCycle_IC(Coords(4,i),:))>1e-6);
            plot(Slopes(RepIndex),LimitCycle_IC(Coords(4,i),RepIndex),'--r','LineWidth',LineWidth);
        end
    else
        h1=plot(Slopes,LimitCycle_IC(Coords(3,:),:),'b','LineWidth',LineWidth);
        hold on
        h2=plot(Slopes,LimitCycle_IC(Coords(4,:),:),'--r','LineWidth',LineWidth);
    end
    if InitCondSubplots==0 || Plots(1)==0
        t1=title('Steady state initial conditions for \theta_{1t}, \theta_{2t} versus slope');
        set(t1,'FontSize',FontSize);
        xh=xlabel(SlopeLabel);
        set(xh,'FontSize',FontSize);
    end
    yh=ylabel('Angular Velocity [rad/s]');
    legend([h1(1);h2(1)],{'\theta_{1t}';'\theta_{2t}'},'Orientation','Horizontal','Location','SouthEast');
    axis([min(Slopes) max(Slopes) -1.25 -0.45]);
    set(gca,'FontSize',FontSize-2);
    set(yh,'FontSize',FontSize);
    if InitCondSubplots==1 && Plots(1)==1
        set(gca,'XTickLabel','');
    end        
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Interactive Phase Plane %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Plots(3)
    X=LimitCycle_X{IPPIndex};
    % Initialize figure
    IPPh_fig=figure(4);
    
    % Set the keyboard input handling function
    set(IPPh_fig,'KeyPressFcn',@KeybHandle);
    
    % Plot
    if InterSubplots==1 && Plots(3)+Plots(5)==2
        subplot(1,2,1);
    end
    IPPh1=plot(X(:,3),X(:,5),'b','LineWidth',LineWidth);
    hold on
    IPPh2=plot(X(:,4),X(:,6),'b','LineWidth',LineWidth);
    axis([-0.5 0.8 -2.8 3.2]);
    xh=xlabel('Angle [rad]');
    yh=ylabel('Angular Velocity [rad/s]');
    IPPt1=title(sprintf('Interactive Phase Plane - Slope=%.2f\\circ',Slopes(IPPIndex)));
    set(gca,'FontSize',FontSize);
    set(xh,'FontSize',FontSize);
    set(yh,'FontSize',FontSize);
    set(IPPt1,'FontSize',FontSize);
end

    %%%%%%%%%%%%%%%%%%%%%%% KeybHandle %%%%%%%%%%%%%%%%%%%%%%%
    % Handle keyboard input
    function KeybHandle(h_obj,evt) %#ok<INUSL>
        if strcmp(evt.Key,'uparrow')
            IPPIndex=IPPIndex+1;
            if IPPIndex>Sections
                IPPIndex=Sections;
            end
        end
        
        if strcmp(evt.Key,'downarrow')
            IPPIndex=IPPIndex-1;
            if IPPIndex<1
                IPPIndex=1;
            end     
        end
        
        if strcmp(evt.Key,'rightarrow')
            IPPIndex=IPPIndex+20;
            if IPPIndex>Sections
                IPPIndex=Sections;
            end
        end
        
        if strcmp(evt.Key,'leftarrow')
            IPPIndex=IPPIndex-20;
            if IPPIndex<1
                IPPIndex=1;
            end     
        end

        X=LimitCycle_X{IPPIndex};
        Slope=Slopes(IPPIndex);
        
        set(IPPh1,'XData',X(:,3));
        set(IPPh1,'YData',X(:,5));
        set(IPPh2,'XData',X(:,4));
        set(IPPh2,'YData',X(:,6));
        set(IPPt1,'String',sprintf('Interactive Phase Plane - Slope=%.2f\\circ',Slope));
        
        if Plots(5)==1 % If interactive PM locus is active, apply change to it as well
            if ~isfield(evt,'Called') % check if this handle was called from the other handle
                evt.Called=1;
                PMKeybHandle([],evt);
            end
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Linearized Poincare Matrix eigenvalues %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Plots(4)
    figure(5);
    subplot(6,1,1:3);
    hold on
    for i=1:5
        plot(Slopes(Region1),abs(LimitCycle_eig(i,Region1)),LineStyles{i},'LineWidth',LineWidth,'Color',Colors{i});
        plot(Slopes(Region2),abs(LimitCycle_eig(i,Region2)),LineStyles{i},'LineWidth',LineWidth,'Color',Colors{i});
        plot(Slopes(Region3),abs(LimitCycle_eig(i,Region3)),LineStyles{i},'LineWidth',LineWidth,'Color',Colors{i});
    end
    
    % Draw bifurcation lines
    X1=Slopes(Region1(end));
    X2=Slopes(Region3(1));
    line([X1 X1],[0 1],'LineWidth',3,'LineStyle','--','Color',[0 0 0]);
    line([X2 X2],[0 1],'LineWidth',3,'LineStyle','--','Color',[0 0 0]);
    
    % Add labels to different areas
    LX1=(Slopes(Region1(end))+Slopes(Region1(1)))/2;
    LX2=(Slopes(Region2(end))+Slopes(Region2(1)))/2;
    LX3=(Slopes(Region3(end))+Slopes(Region3(1)))/2;
    text(LX1,0.96,'I','FontSize',FontSize,'FontName','Times New Roman','FontWeight','bold','HorizontalAlignment','center');
    text(LX2,0.96,'IIa','FontSize',FontSize,'FontName','Times New Roman','FontWeight','bold','HorizontalAlignment','center');
    text(LX3,0.96,'III','FontSize',FontSize,'FontName','Times New Roman','FontWeight','bold','HorizontalAlignment','center');
    
    t1=title('Absolute value of Poincare Map eigenvalues versus slope');
%     xh=xlabel(SlopeLabel);
    yh=ylabel('|\lambda_i|');
    axis([min(Slopes) max(Slopes) -0.05 1]);
    set(gca,'FontSize',FontSize);
%     set(xh,'FontSize',FontSize);
    set(yh,'FontSize',FontSize+2);
    set(t1,'FontSize',FontSize);
    set(gca,'XTickLabel','');
    
    
    subplot(6,1,4:6);
    hold on
    for i=1:5
        plot(Slopes(Region1),abs(LimitCycle_eig(i,Region1)),LineStyles{i},'LineWidth',LineWidth,'Color',Colors{i});
        plot(Slopes(Region2),abs(LimitCycle_eig(i+5,Region2)),LineStyles{i},'LineWidth',LineWidth,'Color',Colors{i});
        plot(Slopes(Region3),abs(LimitCycle_eig(i,Region3)),LineStyles{i},'LineWidth',LineWidth,'Color',Colors{i});
    end
    
    % Draw bifurcation lines
    X1=Slopes(Region1(end));
    X2=Slopes(Region3(1));
    line([X1 X1],[0 1],'LineWidth',3,'LineStyle','--','Color',[0 0 0]);
    line([X2 X2],[0 1],'LineWidth',3,'LineStyle','--','Color',[0 0 0]);
    
    % Add labels to different areas
    LX1=(Slopes(Region1(end))+Slopes(Region1(1)))/2;
    LX2=(Slopes(Region2(end))+Slopes(Region2(1)))/2;
    LX3=(Slopes(Region3(end))+Slopes(Region3(1)))/2;
    text(LX1,0.96,'I','FontSize',FontSize,'FontName','Times New Roman','FontWeight','bold','HorizontalAlignment','center');
    text(LX2,0.96,'IIb','FontSize',FontSize,'FontName','Times New Roman','FontWeight','bold','HorizontalAlignment','center');
    text(LX3,0.96,'III','FontSize',FontSize,'FontName','Times New Roman','FontWeight','bold','HorizontalAlignment','center');
    
%     t1=title('Absolute value of Poincare Map eigenvalues versus slope');
    xh=xlabel(SlopeLabel);
    yh=ylabel('|\lambda_i|');
    axis([min(Slopes) max(Slopes) -0.05 1]);
    set(gca,'FontSize',FontSize);
    set(xh,'FontSize',FontSize);
    set(yh,'FontSize',FontSize+2);
    set(t1,'FontSize',FontSize);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Interactive Poincare Map %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Plots(5)==1
    if InterSubplots==1 && Plots(3)+Plots(5)==2
        figure(IPPh_fig);
        subplot(1,2,2);
    else
        % Initialize figure
        IPMh_fig=figure(6);
        
        % Set the keyboard input handling function
        set(IPMh_fig,'KeyPressFcn',@PMKeybHandle);
    end
        
    % Plot the unit circle
    Circle();
    hold on
    % Select the eigenvalues by index
    switch LimitCycle_Period(IPMIndex)
        case 1
            EigVa=LimitCycle_eig(1:5,IPMIndex);
            EigVb=[];
        case 2
            EigVa=LimitCycle_eig(1:5,IPMIndex);
            EigVb=LimitCycle_eig(6:10,IPMIndex);
        otherwise
            EigVa=[];
            EigVb=[];
    end
    
    % Plot
    IPMha=plot(real(EigVa),imag(EigVa),'*','MarkerSize',10);
    IPMhb=plot(real(EigVa),imag(EigVa),'*r','MarkerSize',10);
    IPMt1=title(sprintf('Numerical PM eigenvalues locus - Slope=%.2f\\circ',Slopes(IPMIndex)));
    xh=xlabel('Real(\lambda)');
    yh=ylabel('Imag(\lambda)');
    IPMtext1=text(0.95,0.95,'Period: 1','FontWeight','bold','HorizontalAlignment','Right');
    IPMaxis=gca;
    
    % Set axis to show the complete unit circle
    axis equal
    if max(abs(LimitCycle_eig(1:5*LimitCycle_Period(IPMIndex),IPMIndex)))<1.05
        set(IPMaxis,'XLim',[-1.05 1.05]);
        set(IPMaxis,'YLim',[-1.05 1.05]);
    else
        set(IPMaxis,'XLimMode','auto');
        set(IPMaxis,'YLimMode','auto');
    end

    set(gca,'FontSize',FontSize);
    set(xh,'FontSize',FontSize);
    set(yh,'FontSize',FontSize);
    set(IPMt1,'FontSize',FontSize);
end

    %%%%%%%%%%%%%%%%%%%%%%% PMKeybHandle %%%%%%%%%%%%%%%%%%%%%%%
    % Handle keyboard input
    function PMKeybHandle(h_obj,evt) %#ok<INUSL>
        if strcmp(evt.Key,'uparrow')
            IPMIndex=IPMIndex+1;
            if IPMIndex>Sections
                IPMIndex=Sections;
            end
        end
        
        if strcmp(evt.Key,'downarrow')
            IPMIndex=IPMIndex-1;
            if IPMIndex<1
                IPMIndex=1;
            end     
        end
        
        if strcmp(evt.Key,'rightarrow')
            IPMIndex=IPMIndex+20;
            if IPMIndex>Sections
                IPMIndex=Sections;
            end
        end
        
        if strcmp(evt.Key,'leftarrow')
            IPMIndex=IPMIndex-20;
            if IPMIndex<1
                IPMIndex=1;
            end     
        end
        
        set(IPMtext1,'string',['Period: ',num2str(LimitCycle_Period(IPMIndex))]);
        % Select the eigenvalues by index
        switch LimitCycle_Period(IPMIndex)
            case 1
                EigVa=LimitCycle_eig(1:5,IPMIndex);
                EigVb=[];
            case 2
                EigVa=LimitCycle_eig(1:5,IPMIndex);
                EigVb=LimitCycle_eig(6:10,IPMIndex);
            otherwise
                EigVa=[];
                EigVb=[];
        end
    
        set(IPMha,'XData',real(EigVa),'YData',imag(EigVa));
        set(IPMhb,'XData',real(EigVb),'YData',imag(EigVb));

        % Set axis to show the complete unit circle
        if max(abs(LimitCycle_eig(1:5*LimitCycle_Period(IPMIndex),IPMIndex)))<1.05
            set(IPMaxis,'XLim',[-1.05 1.05]);
            set(IPMaxis,'YLim',[-1.05 1.05]);
        else
            set(IPMaxis,'XLimMode','auto');
            set(IPMaxis,'YLimMode','auto');
        end

        set(IPMt1,'String',sprintf('Numerical PM eigenvalues locus - Slope=%.2f\\circ',Slopes(IPMIndex)));
        
        if Plots(3)==1 % If interactive phase plane is active, apply change to it as well
            if ~isfield(evt,'Called') % check if this handle was called from the other handle
                evt.Called=1;
                KeybHandle([],evt);
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%% Circle %%%%%%%%%%%%%%%%%%%%%%%
    % Draws a circle
    function [ res ] = Circle()
        Resolution=36;

        coordX=zeros(1,Resolution);
        coordY=zeros(1,Resolution);
        
        for r=1:Resolution
            coordX(1,r)=cos(r/Resolution*2*pi);
            coordY(1,r)=sin(r/Resolution*2*pi);
        end

        h=patch(coordX,coordY,[1 1 1]);
        set(h,'EdgeColor',[0 0 0]);
        set(h,'LineWidth',2,'LineStyle','--');

        res=h;
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Poincare Map eigenvalues locus %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Plots(6)    
    % First we split the slopes into 3 regions
    EigV1=LimitCycle_eig(1:5,Region1)';
    EigV2a=LimitCycle_eig(1:5,Region2)';
    EigV2b=LimitCycle_eig(6:10,Region2)';
    EigV3=LimitCycle_eig(1:5,Region3)';
    
    % We have 3 regions so we'll do 3 subplots
    if LocusSubplots==1
        subplot(4,4,[1 2 5 6]);
    else
        figure(7)
    end
        % Plot the unit circle
        Circle();
        hold on
        line([-1.2 1.2],[0 0],'Color',[0 0 0]);
        line([0 0],[-1.2 1.2],'Color',[0 0 0]);
        % Plot the eigenvalues locus
        for i=1:5
            plot(real(EigV1(:,i)),imag(EigV1(:,i)),'-','LineWidth',LineWidth,'MarkerSize',8,'Color',Colors{i});
        end
        t1=title(sprintf('Region I: %.2f\\circ to %.2f\\circ',Slopes(1),Slopes(Split1-1)));
        yh=ylabel('Imag(\lambda)');
        if LocusSubplots==0
            xh=xlabel('Real(\lambda)');
            set(xh,'FontSize',FontSize);
        end 

        if LocusDirection==1
            b=find(imag(EigV1)==max(max(imag(EigV1))));
            col=floor(b/size(EigV1,1));
            Start=EigV1(b-col*size(EigV1,1),:);
            for jl=1:5
                plot(real(Start(jl)),imag(Start(jl)),'.','MarkerSize',26,'Color',Colors{jl});
            end
            AddDirection(EigV1,3.5);
        end

        % Set axis to show the complete unit circle
        axis([-1.05 1.05 -1.05 1.05]);
        axis image
        axis off
        set(gcf,'Color',[1 1 1]);

        set(gca,'FontSize',FontSize);
        set(yh,'FontSize',FontSize);
        set(t1,'FontSize',FontSize);
    
    if LocusSubplots==1 
        subplot(4,4,[9 10 13 14]);
    else
        figure(8)
    end
        % Plot the unit circle
        Circle();
        hold on
        line([-1.2 1.2],[0 0],'Color',[0 0 0]);
        line([0 0],[-1.2 1.2],'Color',[0 0 0]);
        % Plot the eigenvalues locus
        for i=1:5
            plot(real(EigV2a(:,i)),imag(EigV2a(:,i)),'-','LineWidth',LineWidth,'MarkerSize',8,'Color',Colors{i});
        end
        t1=title(sprintf('Region IIa: %.2f\\circ to %.2f\\circ',Slopes(Split1),Slopes(Split2)));
        xh=xlabel('Real(\lambda)');
        if LocusSubplots==0
            yh=ylabel('Imag(\lambda)');
            set(yh,'FontSize',FontSize);
        end 

        if LocusDirection==1
            Start=EigV2a(1,:);
            for jl=1:5
                plot(real(Start(jl)),imag(Start(jl)),'.','MarkerSize',26,'Color',Colors{jl});
            end
            AddDirection(EigV2a,3);
        end
        
        % Set axis to show the complete unit circle
        axis([-1.05 1.05 -1.05 1.05]);
        axis image
        axis off
        set(gcf,'Color',[1 1 1]);

        set(gca,'FontSize',FontSize);
        set(xh,'FontSize',FontSize);
        set(t1,'FontSize',FontSize);
        
    if LocusSubplots==1 
        subplot(4,4,[11 12 15 16]);
    else
        figure(9)
    end
        % Plot the unit circle
        Circle();
        hold on
        line([-1.2 1.2],[0 0],'Color',[0 0 0]);
        line([0 0],[-1.2 1.2],'Color',[0 0 0]);
        % Plot the eigenvalues locus
        for i=1:5
            plot(real(EigV2b(:,i)),imag(EigV2b(:,i)),'-','LineWidth',LineWidth,'MarkerSize',8,'Color',Colors{i});
        end
        t1=title(sprintf('Region IIb: %.2f\\circ to %.2f\\circ',Slopes(Split1),Slopes(Split2)));
        xh=xlabel('Real(\lambda)');
        if LocusSubplots==0
            yh=ylabel('Imag(\lambda)');
            set(yh,'FontSize',FontSize);
        end 

        if LocusDirection==1
            Start=EigV2b(1,:);
            for jl=1:5
                plot(real(Start(jl)),imag(Start(jl)),'.','MarkerSize',26,'Color',Colors{jl});
            end
            AddDirection(EigV2b,4.5);
        end
        
        % Set axis to show the complete unit circle
        axis([-1.05 1.05 -1.05 1.05]);
        axis image
        axis off
        set(gcf,'Color',[1 1 1]);

        set(gca,'FontSize',FontSize);
        set(xh,'FontSize',FontSize);
        set(t1,'FontSize',FontSize);
    
    if LocusSubplots==1
        subplot(4,4,[3 4 7 8]);
    else
        figure(10)
    end
        % Plot the unit circle
        Circle();
        hold on
        line([-1.2 1.2],[0 0],'Color',[0 0 0]);
        line([0 0],[-1.2 1.2],'Color',[0 0 0]);
        % Plot the eigenvalues locus
        for i=1:5
            plot(real(EigV3(:,i)),imag(EigV3(:,i)),'-','LineWidth',LineWidth,'MarkerSize',8,'Color',Colors{i});
        end
        t1=title(sprintf('Region III: %.2f\\circ to %.2f\\circ',Slopes(Split2+1),Slopes(end)));
        if LocusSubplots==0
            xh=xlabel('Real(\lambda)');
            set(xh,'FontSize',FontSize);
            yh=ylabel('Imag(\lambda)');
            set(yh,'FontSize',FontSize);
        end 

        if LocusDirection==1
            Start=EigV3(1,:);
            for jl=1:5
                plot(real(Start(jl)),imag(Start(jl)),'.','MarkerSize',26,'Color',Colors{jl});
            end
            AddDirection(EigV3,3.5);
        end
        
        % Set axis to show the complete unit circle
        axis([-1.05 1.05 -1.05 1.05]);
        axis image
        axis off
        set(gcf,'Color',[1 1 1]);

        set(gca,'FontSize',FontSize);
        set(t1,'FontSize',FontSize);
    
    
end

    function [] = AddDirection(Values,nArrows)
        [N,M]=size(Values);

        % Find length of locus to evenly distribute arrows
        LocusDist=zeros(1,M);
        for j=1:M
            for a=1:N-1
                LocusDist(j)=LocusDist(j)+sqrt((real(Values(a+1,j))-real(Values(a,j)))^2+(imag(Values(a+1,j))-imag(Values(a,j)))^2);
            end
        end    
        Step=LocusDist/(nArrows+1);
        
        LocusDist=zeros(1,M);
        
        for j=1:M
            for a=1:N-1
                LocusDist(j)=LocusDist(j)+sqrt((real(Values(a+1,j))-real(Values(a,j)))^2+(imag(Values(a+1,j))-imag(Values(a,j)))^2);
                if LocusDist(j)>Step
                    LocusDist(j)=0;
                    
                    Vec1=[real(Values(a,j)) imag(Values(a,j))];
                    Vec2=[real(Values(a+1,j)) imag(Values(a+1,j))];
                    Dir=Vec2-Vec1;
                    Dir=Dir/norm(Dir);
                    Length=0.075;

                    % Plot arrow
                    TrDir=[-Dir(2) Dir(1)];
                    Points=[Vec2; Vec2-Length*Dir-Length/2*TrDir; Vec2-Length*Dir+Length/2*TrDir];
                    patch(Points(:,1), Points(:,2),Colors{j},'EdgeColor',Colors{j});
                    Point4=Vec2-Length*Dir+2*Length*TrDir;
                    if LocusDirNumber==1
                        text(Point4(1),Point4(2),num2str(a),'FontSize',14);
                    end
                end
            end
        end        
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Limit Cycle plots for paper %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Plots(7)==1
    figure(11)
    PlotSlopes=[-10 -5 0 3 7 10];
    SubPlots=[4 7 10 13 16; 5 8 11 14 17; 6 9 12 15 18; 22 25 28 31 34; 23 26 29 32 35; 24 27 30 33 36];
    for i=1:6
        if i==4
            subplot(12,3,SubPlots(i,:));
            CurSlope=find(Slopes>=PlotSlopes(i),1,'first');
            ID=find(abs(diff(LimitCycle_X{CurSlope}(:,6)))>0.1,1,'first');
            plot(LimitCycle_X{CurSlope}(1:ID+1,3),LimitCycle_X{CurSlope}(1:ID+1,5),'--','LineWidth',LineWidth);
            hold on
            plot(LimitCycle_X{CurSlope}(1:ID+1,4),LimitCycle_X{CurSlope}(1:ID+1,6),'r','LineWidth',LineWidth);
            plot(LimitCycle_X{CurSlope}(ID+1:end,3),LimitCycle_X{CurSlope}(ID+1:end,5),'LineWidth',LineWidth);
            plot(LimitCycle_X{CurSlope}(ID+1:end,4),LimitCycle_X{CurSlope}(ID+1:end,6),'--r','LineWidth',LineWidth);
            plot(LimitCycle_X{CurSlope}(1,3),LimitCycle_X{CurSlope}(1,5),'.','MarkerSize',25);
            plot(LimitCycle_X{CurSlope}(1,3),LimitCycle_X{CurSlope}(1,5),'.w','MarkerSize',10);
            plot(LimitCycle_X{CurSlope}(1,4),LimitCycle_X{CurSlope}(1,6),'.r','MarkerSize',25);
            plot(LimitCycle_X{CurSlope}(ID+1,3),LimitCycle_X{CurSlope}(ID+1,5),'.','MarkerSize',25);
            plot(LimitCycle_X{CurSlope}(ID+1,4),LimitCycle_X{CurSlope}(ID+1,6),'.r','MarkerSize',25);
            plot(LimitCycle_X{CurSlope}(ID+1,4),LimitCycle_X{CurSlope}(ID+1,6),'.w','MarkerSize',10);
        else
            subplot(12,3,SubPlots(i,:));
            CurSlope=find(Slopes>=PlotSlopes(i),1,'first');
            plot(LimitCycle_X{CurSlope}(:,3),LimitCycle_X{CurSlope}(:,5),'--','LineWidth',LineWidth);
            hold on
            plot(LimitCycle_X{CurSlope}(:,4),LimitCycle_X{CurSlope}(:,6),'LineWidth',LineWidth);
            plot(LimitCycle_X{CurSlope}(1,3),LimitCycle_X{CurSlope}(1,5),'.','MarkerSize',25);
            plot(LimitCycle_X{CurSlope}(1,3),LimitCycle_X{CurSlope}(1,5),'.w','MarkerSize',10);
            plot(LimitCycle_X{CurSlope}(1,4),LimitCycle_X{CurSlope}(1,6),'.','MarkerSize',25);
        end
        
%         switch i
%             case 4
%                 AddArrows([LimitCycle_X{CurSlope}(:,3), LimitCycle_X{CurSlope}(:,5)], 2, [0 0 1], 5.537, 0);
%             case 5
%                 AddArrows([LimitCycle_X{CurSlope}(:,3), LimitCycle_X{CurSlope}(:,5)], 6, [0 0 1], 5.537, [1 2 3 5]);
%             case 6
%                 AddArrows([LimitCycle_X{CurSlope}(:,3), LimitCycle_X{CurSlope}(:,5)], 6, [0 0 1], 5.537, [1 2 3 5]);
%             otherwise
%                 AddArrows([LimitCycle_X{CurSlope}(:,3), LimitCycle_X{CurSlope}(:,5)], 5, [0 0 1], 5.537, [2 3 5]);
%         end
        
        axis([-0.5 0.62 -2.8 3.4]);
        set(gca,'FontSize',FontSize);
        if Slopes(CurSlope)<=0
            title(['Slope: ',num2str(Slopes(CurSlope)),'\circ'],'FontSize',FontSize);
            set(gca,'XTickLabel','');
        else
            title(['Slope: +',num2str(Slopes(CurSlope)),'\circ'],'FontSize',FontSize);
        end
        if i~=1 && i~=4
            set(gca,'YTickLabel','');
        end            
    end
    
    subplot(12,3,SubPlots(5,:));
    xlabel('Leg angle [rad]','FontSize',FontSize);
%     subplot(2,3,SubPlots(1,:));
%     ylabel('Leg angular velocity [rad/s]                                                                                ','FontSize',FontSize);
    
    hax = axes('Position', [0, 0, 1, 1],'Visible','off');
    text(0.05,0.5,'Leg angular velocity [rad/s]','FontSize',FontSize,'Units','normalized','rotation',90,'HorizontalAlignment','center');
    text(0.5,0.99,'Limit cycles generated by the Full Feedback Controller for different slopes',...
        'FontSize',FontSize,'Units','normalized','HorizontalAlignment','center','VerticalAlignment','top');
end

%     function [] = AddArrows(Values,nArrows,Color,AR,skip)
%         N=length(Values);
% 
%         % Find length of locus to evenly distribute arrows
%         LocusDist=0;
%         for lo=1:N-1
%             LocusDist=LocusDist+sqrt((Values(lo+1,1)-Values(lo,1))^2+(Values(lo+1,2)-Values(lo,2))^2);
%         end
%         Step=LocusDist/(nArrows+1);
%         
%         LocusDist=0;
%         count=0;        
%         for lo=1:N-1
%             LocusDist=LocusDist+sqrt((Values(lo+1,1)-Values(lo,1))^2+(Values(lo+1,2)-Values(lo,2))^2);
%             if LocusDist>Step
%                 LocusDist=0;
%                 
%                 count=count+1;
%                 if ~isempty(find(count==skip,1,'first'))
%                     continue
%                 end
% 
%                 Vec1=[Values(lo,1) Values(lo,2)];
%                 Vec2=[Values(lo+1,1) Values(lo+1,2)];
%                 Dir=Vec2-Vec1;
%                 Dir=Dir/norm(Dir);
%                 Length=0.1*Step;
% 
%                 % Plot arrow
%                 TrDir=[-Dir(2) Dir(1)];
%                 Points=[Vec2; Vec2-Length*Dir-Length/2*TrDir; Vec2-Length*Dir+Length/2*TrDir];
%                 
%                 % Correct aspect ratio
%                 CoM=mean(Points,1);
%                 Points=Points-repmat(CoM,3,1);
%                 Points(:,2)=Points(:,2)*AR;
%                 Points=Points+repmat(CoM,3,1);
%                 
%                 patch(Points(:,1), Points(:,2),Color,'EdgeColor',Color);
%             end
%         end   
%     end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Compass Gait Model Figure %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Plots(8)==1
    Robot=CompassBiped();
    InitCond=[1.88*pi, 0.2*pi, 0, 0];
    Robot.CircRes=72;
    Robot.LinkRes=24;
    Robot.RenderParams=1;
    
    Floor=Terrain(0,1,1);
    figure(12)
    Robot.Render(InitCond);
    Floor.Render(-0.4,1.2);
    
    axis([-0.4 1.2 -0.5 1.5])
    axis off
    set(gcf,'Color',[1 1 1])
    
    set(0,'Units','Pixels');
    ScSize=get(0,'ScreenSize');
    set(12,'Units','Pixels');
    set(12,'OuterPosition',ScSize+[600 200 -1200 -325]);
    addpath('ExportFig');
    export_fig 'CompassBiped' -r600 -painters -tif -nocrop -transparent
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Interactive Ground Reaction Forces Plot %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Plots(9)==1
    Robot=CompassBiped();
    
    ThisGRF=GRF{IGRFIndex};
    ThisMeasGRF=MeasuredGRF{IGRFIndex};
    ThisT=LimitCycle_T{IGRFIndex};
    
    % Initialize figure
    IGRFh_fig=figure(4);
    
    % Set the keyboard input handling function
    set(IGRFh_fig,'KeyPressFcn',@KeybHandleGRF);
    
    % Plot GRF Fy
    subplot(3,1,1);
    IGRFh1=plot(ThisT,ThisGRF(2,:),'LineWidth',LineWidth);
    hold on
    plot(ThisT,Robot.GetWeight()*ones(size(ThisT)),'--k','LineWidth',LineWidth);
    yh=ylabel('GRF F_y [N]');
    IGRFt1=title(sprintf('Interactive GRF Plot - Slope=%.2f\\circ',Slopes(IGRFIndex)));
    set(gca,'FontSize',FontSize);
    set(yh,'FontSize',FontSize);
    set(IGRFt1,'FontSize',FontSize);
    axis([ThisT(1) ThisT(end) 0.995*min(ThisGRF(2,:)) 1.005*max(ThisGRF(2,:))]);

    % Plot GRF Fx
    subplot(3,1,2);
    IGRFh2=plot(ThisT,ThisGRF(1,:),'LineWidth',LineWidth);
    yh=ylabel('GRF F_x [N]');
    set(gca,'FontSize',FontSize);
    set(yh,'FontSize',FontSize);
    axis([ThisT(1) ThisT(end) 1.05*min(ThisGRF(1,:)) 1.05*max(ThisGRF(1,:))]);
    
    % Plot Measured GRF
    subplot(3,1,3);
    IGRFh3=plot(ThisT,ThisMeasGRF,'LineWidth',LineWidth);
    hold on
    IGRFh4=plot(ThisT,MeanMeasuredGRF(IGRFIndex)*ones(size(ThisT)),'--k','LineWidth',LineWidth);
    xh=xlabel('Time [sec]');
    yh=ylabel('GRF MEasured [N]');
    set(gca,'FontSize',FontSize);
    set(xh,'FontSize',FontSize);
    set(yh,'FontSize',FontSize);
    axis([ThisT(1) ThisT(end) 0.995*min(ThisMeasGRF) 1.005*max(ThisMeasGRF)]);
end

    %%%%%%%%%%%%%%%%%%%%%%% KeybHandle %%%%%%%%%%%%%%%%%%%%%%%
    % Handle keyboard input
    function KeybHandleGRF(h_obj,evt) %#ok<INUSL>
        if strcmp(evt.Key,'uparrow')
            IGRFIndex=IGRFIndex+1;
            if IGRFIndex>Sections
                IGRFIndex=Sections;
            end
        end
        
        if strcmp(evt.Key,'downarrow')
            IGRFIndex=IGRFIndex-1;
            if IGRFIndex<1
                IGRFIndex=1;
            end     
        end
        
        if strcmp(evt.Key,'rightarrow')
            IGRFIndex=IGRFIndex+20;
            if IGRFIndex>Sections
                IGRFIndex=Sections;
            end
        end
        
        if strcmp(evt.Key,'leftarrow')
            IGRFIndex=IGRFIndex-20;
            if IGRFIndex<1
                IGRFIndex=1;
            end     
        end

        ThisGRF=GRF{IGRFIndex};
        ThisMeasGRF=MeasuredGRF{IGRFIndex};
        ThisT=LimitCycle_T{IGRFIndex};
    
        set(IGRFh1,'XData',ThisT);
        set(IGRFh1,'YData',ThisGRF(2,:));
        set(IGRFh2,'XData',ThisT);
        set(IGRFh2,'YData',ThisGRF(1,:));
        set(IGRFh3,'XData',ThisT);
        set(IGRFh3,'YData',ThisMeasGRF);
        set(IGRFh4,'XData',ThisT);
        set(IGRFh4,'YData',MeanMeasuredGRF(IGRFIndex)*ones(size(ThisT)));
        set(IGRFt1,'String',sprintf('Interactive GRF Plot - Slope=%.2f\\circ',Slopes(IGRFIndex)));
        
        subplot(3,1,1);
        axis([ThisT(1) ThisT(end) 0.995*min(ThisGRF(2,:)) 1.005*max(ThisGRF(2,:))]);
        subplot(3,1,2);
        axis([ThisT(1) ThisT(end) 1.05*min(ThisGRF(1,:)) 1.05*max(ThisGRF(1,:))]);
        subplot(3,1,3);
        axis([ThisT(1) ThisT(end) 0.995*min(ThisMeasGRF) 1.005*max(ThisMeasGRF)]);
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Energy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Plots(10)==1
    figure
    plot(Slopes,StepLength);
    hold on
    plot(DPSlopes,DPStepLength,'r');
    title('Step length');
    figure
    plot(Slopes,dKineticEIm);
    hold on
    plot(DPSlopes,DPdKineticEIm,'r');
    title('Kinetic energy lost at impact');
    figure
    plot(Slopes,dKineticEFr);
    hold on
    plot(DPSlopes,DPdKineticEFr,'r');
    title('Kinetic energy lost by friction');
    figure
    plot(Slopes,dPotentialE);
    hold on
    plot(DPSlopes,DPdPotentialE,'r');
    title('Potential energy lost/gained per step');
    figure
    plot(Slopes,ControllerE);
    hold on
    plot(DPSlopes,DPControllerE,'r');
    title('Controller energy input');
    figure
    plot(Slopes,ControllerE-dPotentialE+dKineticEIm+dKineticEFr);
    title('Total energy change');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Mega Interactive Plot %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Plots(11)==1
    Slope = Slopes(IMIPIndex);
    X = LimitCycle_X{IMIPIndex};
    T = LimitCycle_T{IMIPIndex};
    Torques = LimitCycle_Torques{IMIPIndex};
    GRF = LimitCycle_GRF{IMIPIndex};
    ZMP = LimitCycle_ZMP{IMIPIndex};

    % Handles
    IMIPh1 = zeros(2,1);
    IMIPh2 = zeros(2,1);
    IMIPh3 = 0;
    IMIPh4 = zeros(2,1);
    IMIPh5 = zeros(2,1);
    IMIPh6 = 0;
    
    % Initialize figure
    IMIP_fig = figure();
    
    % Set the keyboard input handling function
    set(IMIP_fig,'KeyPressFcn',@KeybHandleMIP);
    
    subplot(3,2,1);
    % Theta 1, Theta 2
    IMIPh1(1) = plot(T,X(:,3));
    hold on
    IMIPh1(2) = plot(T,X(:,4),'r');
    legend('\theta_1','\theta_2');
    IMIPt1 = title(['System state time series - Slope: ',num2str(Slope)]);
    
    subplot(3,2,3);
    % Theta 1 dot, Theta 2 dot
    IMIPh2(1) = plot(T,X(:,5));
    hold on
    IMIPh2(2) = plot(T,X(:,6),'r');
    legend('\theta_1t','\theta_2t');
    
    subplot(3,2,5);
    % CPG phase phi
    IMIPh3 = plot(T,X(:,7));
    legend('\phi');
    
    subplot(3,2,2);
    % Torques
    IMIPh4(1) = plot(T,Torques(:,1));
    hold on
    IMIPh4(2) = plot(T,Torques(:,2),'r');
    legend('Swing','Stance');
    title('Applied torques');
    
    subplot(3,2,4);
    % Ground reaction forces
    IMIPh5(1) = plot(T,GRF(:,1));
    hold on
    IMIPh5(2) = plot(T,GRF(:,2),'r');
    legend('GRF_x','GRF_y');
    title('Ground reaction forces');
    
    subplot(3,2,6);
    % ZMP
    IMIPh6 = plot(T,ZMP);
    legend('ZMP_x');
    title('ZMP position');
end

    %%%%%%%%%%%%%%%%%%%%%%% KeybHandle %%%%%%%%%%%%%%%%%%%%%%%
    % Handle keyboard input
    function KeybHandleMIP(h_obj,evt) %#ok<INUSL>
        if strcmp(evt.Key,'uparrow')
            IMIPIndex=IMIPIndex+1;
            if IMIPIndex>Sections
                IMIPIndex=Sections;
            end
        end
        
        if strcmp(evt.Key,'downarrow')
            IMIPIndex=IMIPIndex-1;
            if IMIPIndex<1
                IMIPIndex=1;
            end     
        end
        
        if strcmp(evt.Key,'rightarrow')
            IMIPIndex=IMIPIndex+20;
            if IMIPIndex>Sections
                IMIPIndex=Sections;
            end
        end
        
        if strcmp(evt.Key,'leftarrow')
            IMIPIndex=IMIPIndex-20;
            if IMIPIndex<1
                IMIPIndex=1;
            end     
        end

        Slope = Slopes(IMIPIndex);
        X = LimitCycle_X{IMIPIndex};
        T = LimitCycle_T{IMIPIndex};
        Torques = LimitCycle_Torques{IMIPIndex};
        GRF = LimitCycle_GRF{IMIPIndex};
        ZMP = LimitCycle_ZMP{IMIPIndex};
    
        subplot(3,2,1);
        UpdatePlot(IMIPh1,T,X(:,3:4));
        subplot(3,2,3);
        UpdatePlot(IMIPh2,T,X(:,5:6));
        subplot(3,2,5);
        UpdatePlot(IMIPh3,T,X(:,7));
        
        subplot(3,2,2);
        UpdatePlot(IMIPh4,T,Torques);
        subplot(3,2,4);
        UpdatePlot(IMIPh5,T,GRF);
        subplot(3,2,4);
        UpdatePlot(IMIPh6,T,ZMP);
        
        set(IMIPt1,'String',['System state time series - Slope: ',num2str(Slope)]);      
    end

    function UpdatePlot(h,T,X)
        for d=1:length(h)
            set(h(d),'XData',T);
            set(h(d),'YData',X(:,d));
        end
        
        MinT = min(T);
        MaxT = max(T);
        MinX = min(min(X));
        MaxX = max(max(X));
        margin = 0.05*abs(MaxX-MinX);
        MinX = MinX-margin;
        MaxX = MaxX+margin;
        axis([MinT, MaxT, MinX, MaxX]);
    end
     

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Mega Interactive Plot %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Plots(12)==1
    figure
    plot(Slopes,MaxFriction);
    xlabel('Slopes [deg]');
    ylabel('Required friction coefficient');
    title('Required friction VS slope');
    
    figure
    plot(Slopes,MaxZMP);
    xlabel('Slopes [deg]');
    ylabel('ZMP position');
    title('Required foot length');
end

%% Closed Loop Controller Parameters

if Plots(13)==1
    figure(5);
    % Closed Loop - Muscle Torques
    plot(SlopesCL,TorquesCL,'*');
    title('Closed Loop - Muscle Torques');
    xlabel(SlopeLabel);
    ylabel('Torque [Nm]');
    legend('Extensor','Flexor','Location','SouthWest');
    axis([min(SlopesCL) max(SlopesCL) -7 3]);
end


%% Potential, Kinetic and Total Energy

if Plots(14)==1
    % Potential, Kinetic and Total Energy
    EnergiesCL=zeros(Sections(CL),4*4);
    for i=1:Sections(CL)
        EnergiesCL(i,1:4)=GetEnergies(IC_SS_CL{i}.InitCond(1,:),SlopesCL(i));
        EnergiesCL(i,5:8)=GetEnergies(IC_SS_CL{i}.InitCond(2,:),SlopesCL(i));
        EnergiesCL(i,9:12)=GetEnergies(IC_SS_CL{i}.InitCond(3,:),SlopesCL(i));
        EnergiesCL(i,13:16)=GetEnergies(IC_SS_CL{i}.InitCond(4,:),SlopesCL(i));
    end
    
    figure(7);
    h1=plot(SlopesCL,EnergiesCL(:,[1 5 9 13]),'*b');
    hold on
    h2=plot(SlopesCL,EnergiesCL(:,[2 6 10 14]),'*r');
    h3=plot(SlopesCL,EnergiesCL(:,[3 7 11 15]),'*k');
    h4=plot(SlopesCL,EnergiesCL(:,[4 8 12 16]),'sb');
    title('Closed Loop - KE, PE, \DeltaPE and Total Energy');
    xlabel(SlopeLabel);
    ylabel('Energy [Joule]');
    legend([h1(1);h2(1);h3(1);h4(1)],{'KE';'PE';'\DeltaPE';'TotalE'});
end

%% Off Center of Mass distance

if Plots(15)
    % Closed Loop - Off COM Distance
    figure(9);
    plot(SlopesCL,OCOMDCL,'*b');
    title('Closed Loop - Off COM Distance');
    xlabel(SlopeLabel);
    ylabel('Distance [%]');
    axis([min(SlopesCL) max(SlopesCL) -60 100]);
end


%% Off COM distance by leg apperture

if Plots(16)
    alpha=linspace(3*pi/180,45*pi/180,20);
    OffCenter1=zeros(size(alpha));
    Slope=2;
    for i=1:length(alpha)
        X=[-alpha(i)/2+Slope*RMult,alpha(i)/2+Slope*RMult];
        OffCenter1(i)=OffCOM(X,Slope);
    end
    OffCenter2=zeros(size(alpha));
    Slope=4;
    for i=1:length(alpha)
        X=[-alpha(i)/2+Slope*RMult,alpha(i)/2+Slope*RMult];
        OffCenter2(i)=OffCOM(X,Slope);
    end
    OffCenter3=zeros(size(alpha));
    Slope=6;
    for i=1:length(alpha)
        X=[-alpha(i)/2+Slope*RMult,alpha(i)/2+Slope*RMult];
        OffCenter3(i)=OffCOM(X,Slope);
    end
    OffCenter4=zeros(size(alpha));
    Slope=8;
    for i=1:length(alpha)
        X=[-alpha(i)/2+Slope*RMult,alpha(i)/2+Slope*RMult];
        OffCenter4(i)=OffCOM(X,Slope);
    end
    OffCenter5=zeros(size(alpha));
    Slope=10;
    for i=1:length(alpha)
        X=[-alpha(i)/2+Slope*RMult,alpha(i)/2+Slope*RMult];
        OffCenter5(i)=OffCOM(X,Slope);
    end
    OffCenter6=zeros(size(alpha));
    Slope=12;
    for i=1:length(alpha)
        X=[-alpha(i)/2+Slope*RMult,alpha(i)/2+Slope*RMult];
        OffCenter6(i)=OffCOM(X,Slope);
    end

    figure(10);
    plot(alpha,OffCenter1);
    hold on
    plot(alpha,OffCenter2,'--');
    plot(alpha,OffCenter3,'.-');
    plot(alpha,OffCenter4);
    plot(alpha,OffCenter5,'--');
    plot(alpha,OffCenter6,'.-');
    plot(alpha,100*ones(size(alpha)),'r');
    title('Off COM Distance by Slope and Leg Apperture');
    xlabel('Leg apperture [rad]');
    ylabel('Distance [%]');
    legend('2\circ','4\circ','6\circ','8\circ','10\circ','12\circ');
end


%% Energy Loss by Leg Apperture

% Nondimensional Parameters
Robot = CompassBiped();
% alphaND=mh/m;
% Igag=I/m/L^2;
% M1=5/2+2*alphaND+2*Igag;
% M2=1/2+2*Igag;

    function [ELoss]=EnergyLoss(X)
        EnergyBefore=GetEnergies(X,Slope);
        
        cSNS=cos(X(1)-X(2));
        % Calculate Impact
        P=[            1-4*Igag             0     ;
            -(3+4*alphaND)*cSNS-4*Igag  1-4*Igag  ];

        Q=[    2*cSNS       -2*M2    ;
            2*cSNS-2*M1  2*cSNS-2*M2 ];
        Xdotnew=pinv(Q)*P*[X(3); X(4)];
        X(3)=Xdotnew(1);
        X(4)=Xdotnew(2);
        EnergyAfter=GetEnergies(X,Slope);
        
        ELoss=EnergyBefore(1)-EnergyAfter(1);
    end

if Plots(17)
    alpha=linspace(3*pi/180,45*pi/180,20);
    ELoss1=zeros(size(alpha));
    Slope=0;
    for i=1:length(alpha)
        X=[-alpha(i)/2+Slope*RMult,alpha(i)/2+Slope*RMult,-0.8,-0.6];
        ELoss1(i)=EnergyLoss(X);
    end

    figure(15);
    plot(alpha,ELoss1);
%     hold on
%     plot(alpha,OffCenter2,'--');
%     plot(alpha,OffCenter3,'.-');
%     plot(alpha,OffCenter4);
%     plot(alpha,OffCenter5,'--');
%     plot(alpha,OffCenter6,'.-');
%     plot(alpha,100*ones(size(alpha)),'r');
    title('Energy Loss by Leg Apperture');
    xlabel('Leg apperture [rad]');
    ylabel('Energy Lost [N]');
%     legend('2\circ','4\circ','6\circ','8\circ','10\circ','12\circ');
end


%% Phase Planes for paper

if Plots(18)
    figure(12);
    
    %       [left bottom ?? top]
%     myinset=[10000000 10 0 10];
    myinset1=[0 0.1 0.33 0.9];
    myinset2=[0.36 0.1 0.66 0.9];
    myinset3=[0.69 0.1 1 0.9];
    
    X=XTempOL{1};
    subplot(1,3,1);
    h1=plot(X(:,1),X(:,3),'r','LineWidth',2);
    hold on
    h2=plot(X(:,2),X(:,4),'r','LineWidth',2);
    axis([-0.5 0.8 -2.4 2.8]);
%     set(gca,'LooseInset',myinset);
    set(gca,'Position',myinset1)
    yh=ylabel('Angular Velocity [rad/s]');
    t1=title(sprintf('Slope=%.2f\circ',SlopesOL(1)));
    set(gca,'FontSize',14);
    
    X=XTempOL{5};
    subplot(1,3,2);
    h1=plot(X(:,1),X(:,3),'r','LineWidth',2);
    hold on
    h2=plot(X(:,2),X(:,4),'r','LineWidth',2);
    axis([-0.5 0.8 -2.4 2.8]);
%     set(gca,'LooseInset',myinset);
    set(gca,'Position',myinset2)
    xh=xlabel('Angle [rad]');
    t2=title(sprintf('Slope=%.2f\circ',SlopesOL(5)));
    set(gca,'FontSize',14);
    
    X=XTempOL{29};
    subplot(1,3,3);
    h1=plot(X(:,1),X(:,3),'r','LineWidth',2);
    hold on
    h2=plot(X(:,2),X(:,4),'r','LineWidth',2);
    axis([-0.5 0.8 -2.4 2.8]);
%     set(gca,'LooseInset',myinset);
    set(gca,'Position',myinset3)
    t3=title(sprintf('Slope=%.2f\circ',SlopesOL(29)));
    set(gca,'FontSize',14);
    
    set(xh,'FontSize',14);
    set(yh,'FontSize',14);
    set(t1,'FontSize',14);
    set(t2,'FontSize',14);
    set(t3,'FontSize',14);
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Stick-figure images %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Plots(19)==1
    [T1,X1]=Simulation( 3, 0, 0.7, 0, 2 );
    [T2,X2]=Simulation( 3, 0, -0.9, 0, 2 );
    [T3,X3]=Simulation( 3, 0, 5.5, 1, 2 );
    [T4,X4]=Simulation( 3, 0, -4, 1, 2 );
    [T5,X5]=Simulation( 3, 0, 9.5, 2, 2 );
    [T6,X6]=Simulation( 3, 0, -9.5, 2, 2 );

    MaxX=max(X5(end,1),X6(end,1));
    MaxY=X5(end,2);
    MinY=X6(end,2);
    StickFigures( T1, X1, 0, 0.7, 0, 'OL_Up.png', MaxX, MaxY );
    StickFigures( T2, X2, 0, -0.9, 0, 'OL_Down.png', MaxX, MinY );
    StickFigures( T3, X3, 0, 5.5, 1, 'Partial_Up.png', MaxX, MaxY );
    StickFigures( T4, X4, 0, -4, 1, 'Partial_Down.png', MaxX, MinY );
    StickFigures( T5, X5, 0, 9.5, 2, 'Full_Up.png', MaxX, MaxY );
    StickFigures( T6, X6, 0, -9.5, 2, 'Full_Down.png', MaxX, MinY );
end

%% Auxiliary Functions

    %%%%%%%%%%%%%%%%%%%%%%% Adaptation %%%%%%%%%%%%%%%%%%%%%%%
    % Neural adaptation function
    function [ETorque,FTorque]=Adaptation(Slope)        
        Ek=2/m/g/L;
        
        Phi=Slope*RMult;
               
        if Phi>0
            omega=omega0+omega0/(2*pi)*2.25*Phi;
            EEnergy=EEnergy0+7*Phi;
            FEnergy=FEnergy0-150*Phi;
        else
            omega=omega0+omega0/(2*pi)*5*Phi;%+150*Phi^3;
            EEnergy=EEnergy0+7*Phi;
            FEnergy=FEnergy0-120*Phi;
        end
        
        delta_t1=0.0995861/omega;
        delta_t2=0.3983445/omega;
    
        ETorque=m*L^2/2*Ek*EEnergy/delta_t1;
        FTorque=m*L^2/2*Ek*FEnergy/delta_t2;
    end

    %%%%%%%%%%%%%%%%%%%%%%% OffCOM %%%%%%%%%%%%%%%%%%%%%%%
    % Calculates the distance from the Robot's center of
    % mass (COM) to the center of the support base in %
    function [Distance]=OffCOM(X,Slope)
        % Calculate COM height distance from floor (perpendicular)
        [COMx, COMy]=RobotCOM(X-Slope*RMult);
        Distance=COMy*tan(Slope*RMult)/abs(COMx)*100;
        
%         % Method 2
%         % Calculate the COM's projection on the ground
%         COMProj=[0;0];
%         COMProj(1)=RobotCOM(X);
%         COMProj(2)=tan(Slope*RMult)*COMProj(1);
%         
%         % Direction of movement
%         [x,y]=GetPos(X,4); % Non-support foot
%         FootDir=[x;y];
%         MoveDir=[1;tan(Slope*RMult)];
%         Factor=sign(FootDir'*MoveDir);
%         Middle=FootDir/2;
%         Dist=COMProj-Middle;
%         if Dist*FootDir'>0
%             Distance=-Factor*norm(Dist)/norm(Middle)*100;
%         else
%             Distance=norm(Dist)/norm(Middle)*100;
%         end
        
%         % Method 1
%         HalfStepLength=norm(FootDir)/2;
%         if FootDir*COMProj'<0
%             % COM falls outside of support leg
%             Distance=-1-norm(COMProj)/HalfStepLength;
%         else
%             Distance=abs(norm(COMProj)-HalfStepLength)/HalfStepLength;
%         end
%         Distance=Distance*100;
    end

    %%%%%%%%%%%%%%%%%%%%%%% RobotCOM %%%%%%%%%%%%%%%%%%%%%%%
    % Calculates the Robot's center of mass (COM)
    % relative to the support foot
    function [COMx,COMy]=RobotCOM(X)
        [M1x,M1y]=GetPos(X,1); % Non-support mass
        [M2x,M2y]=GetPos(X,2); % Support mass
        [M3x,M3y]=GetPos(X,3); % Hip mass
        
        COMx=(m*M1x+m*M2x+mh*M3x)/(m+m+mh);
        COMy=(m*M1y+m*M2y+mh*M3y)/(m+m+mh);
    end

    %%%%%%%%%%%%%%%%%%%%%%% GetPos %%%%%%%%%%%%%%%%%%%%%%%
    % Returns the position of each part of the robot
    function [x, y] = GetPos(X,which)
        sS=sin(X(1)); cS=cos(X(1));
        sNS=sin(X(2)); cNS=cos(X(2));
        
        if which==1 % Non-support mass
            x=-L*sS+L/2*sNS;
            y=L*cS-L/2*cNS;
        end
        if which==2 % Support mass
            x=-L/2*sS;
            y=L/2*cS;
        end
        if which==3 % Hip mass
            x=-L*sS;
            y=L*cS;
        end
        if which==4 % Non-support foot
            x=-L*sS+L*sNS;
            y=L*cS-L*cNS;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%% GetVel %%%%%%%%%%%%%%%%%%%%%%%
    % Returns the velocity of each part of the robot
    function [x, y] = GetVel(X,which)
        sS=sin(X(1)); cS=cos(X(1));
        sNS=sin(X(2)); cNS=cos(X(2));
        qSt=X(3);
        qNSt=X(4);
        
        if which==1 % Non-support mass
            x=-L*cS*qSt*+L/2*cNS*qNSt;
            y=-L*sS*qSt+L/2*sNS*qNSt;
        end
        if which==2 % Support mass
            x=-L/2*cS*qSt;
            y=-L/2*sS*qSt;
        end
        if which==3 % Hip mass
            x=-L*cS*qSt;
            y=-L*sS*qSt;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%% GetEnergies %%%%%%%%%%%%%%%%%%%%%%%
    % Returns the different energies of the robot
    function [SysEnergy] = GetEnergies(X,Slope)
        SysEnergy=zeros(1,4); % KE, PE, DeltaPE & TotalE
        qSt=X(3);
        qNSt=X(4);
        
        [M1x,M1y]=GetPos(X,1); % Non-support mass position
        [M2x,M2y]=GetPos(X,2); % Support mass position
        [M3x,M3y]=GetPos(X,3); % Hip mass position
        
        [M1vx,M1vy]=GetVel(X,1); % Non-support mass velocity
        [M2vx,M2vy]=GetVel(X,2); % Support mass velocity
        [M3vx,M3vy]=GetVel(X,3); % Hip mass velocity
        
        SysEnergy(1)=1/2*m*(M1vx^2+M1vy^2)+1/2*m*(M2vx^2+M2vy^2)+1/2*mh*(M3vx^2+M3vy^2)+1/2*I*(qSt^2+qNSt^2);
        SysEnergy(2)=m*g*(M1y+M2y)+mh*g*M3y;
        
        StepLength=2*L*sin(abs(X(1)-X(2))/2);
        SysEnergy(3)=(m+m+mh)*g*StepLength*sin(Slope*RMult);
        
        % In order to get the Energy Input from the motor
        % we must calculate Int(Torque*theta_dot,t=t0..t1)
        
        % In order to get the Energy Loss from impact
        % we must calculate the velocities prior to impact
        % and obtain DeltaKE
        
        SysEnergy(4)=sum(SysEnergy(1:3));
    end

end
