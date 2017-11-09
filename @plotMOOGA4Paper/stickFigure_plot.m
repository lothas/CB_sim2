function [Sim] = stickFigure_plot(obj,whichGA,geneNum, GenID, Dur, timestep, filename)
%PLOT_SEQ plots the CPG output for a given sequence

GA = obj.data{1,whichGA}.GA;

Generation = geneNum;

% % Animation properties
% FPS = 25;
% Background = [0 0 0]; % black
% % Background = [1 1 1]; % white
% 
% % Render handles
% AnFig = 0; tCOM = 0; COMx0 = 0; COMy0 = 0;
% AnWinW = 0.4; AnWinH = 0.8; AR = 0;
% FlMin = 0; FlMax = 0; HeightMin = 0; HeightMax = 0;

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

Sim.Graphics = 0;
Sim.Follow = 1;

Sim = GA.Gen.Decode(Sim, GA.Seqs(GenID,:,Generation));

% disp(GA.Gen.seq2str(GA.Seqs(GenID,:,Generation)));

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

%Display X
figure(1);
hold on;
x=0;
cnt1=0;
N=2;
supportPos = diff(Sim.Out.SuppPos(:,1));
alpha = linspace(1,0.2,(length(Sim.Out.X)-1)); %0.2
leftLeg = [0,0,1];
rightLeg = [1,0,0];
 
for i=1:(length(Sim.Out.X)-1)
    if supportPos(i,1)
        x = Sim.Out.SuppPos(i,1);
        leftLeg = abs([1,0,1]-leftLeg);
        rightLeg = abs([1,0,1]-rightLeg);
        DrawBiPedRobot2(Sim.Out.X(i,:),leftLeg,rightLeg,x,alpha(1,i))
    else
        cnt1=cnt1+1;
        if(cnt1>=N)
              DrawBiPedRobot2(Sim.Out.X(i,:),leftLeg,rightLeg,x,alpha(1,i)) 
              cnt1=0;
        end
    end
   
    
end

end

function [] = DrawBiPedRobot2(X,leftLeg,rightLeg,x,alpha)
 
 figure(gcf);
 hold on;
 
 L=1;
 b=0.01;
 
 X1 = [0+b*cos(X(1)),0-b*cos(X(1)),0-L*sin(X(1))-b*cos((X(1))),0-L*sin(X(1))+b*cos((X(1)))]+x;
 Y1 = [0+b*sin(X(1)),0-b*sin(X(1)),0+L*cos(X(1))-b*sin((X(1))),0+L*cos(X(1))+b*sin((X(1)))];
 fill(X1,Y1,leftLeg,'EdgeColor','k','FaceAlpha',alpha)
% patch(X1,Y1,[0 0 0 0],[0 0 1]);
 
 X2 = [-L*sin(X(1))-b*cos(X(2)),-L*sin(X(1))+b*cos(X(2)),-L*sin(X(1))+L*sin(X(2))+b*cos((X(2))),-L*sin(X(1))+L*sin(X(2))-b*cos((X(2)))]+x;
 Y2 = [L*cos(X(1))-b*sin(X(2)),L*cos(X(1))+b*sin(X(2)),L*cos(X(1))-L*cos(X(2))+b*sin((X(2))),L*cos(X(1))-L*cos(X(2))-b*sin((X(2)))];
 fill(X2,Y2,rightLeg,'EdgeColor','k','FaceAlpha',alpha)
%  axis equal


end