function [Sim] = CBstick_Figure_plot(obj,whichGA,geneNum, GenID, Dur, timestep, filename)
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

% Load Variables:
X = Sim.Out.X;


%Display X
figure(1);
hold on;

supportPos = diff(Sim.Out.SuppPos(:,1));
alpha = linspace(1,0.2,(length(Sim.Out.X)-1)); %0.2
leftLeg = [0,0,1];
rightLeg = [1,0,0];
 

for i=1:(length(X)-1)
    
    if i==1 || ~mod(i,10)
        xS = Sim.Out.SuppPos(i,1);
        yS = Sim.Out.SuppPos(i,2);

        [x1,y1] = GetLegsPos(Sim.Mod,X(i,:),xS,yS,'Hip');
        [x2,y2] = GetLegsPos(Sim.Mod,X(i,:),xS,yS,'NS');
        [x1a,y1a] = GetLegsPos(Sim.Mod,X(i,:),xS,yS,'m1');
        [x2a,y2a] = GetLegsPos(Sim.Mod,X(i,:),xS,yS,'m2');

        % Render masses as circles
        DrawCircle(Sim.Mod, x1a, y1a, 0, Sim.Mod.m_radius, Sim.Mod.m_color,1);
        DrawCircle(Sim.Mod, x2a, y2a, 0, Sim.Mod.m_radius, Sim.Mod.m_color,2);
        DrawCircle(Sim.Mod, x1, y1, 0, Sim.Mod.mh_radius, Sim.Mod.mh_color,3);
        % Render links
        DrawLink(Sim.Mod, xS, yS, x1, y1, 0);
        DrawLink(Sim.Mod, x1, y1, x2, y2, 0);

        if i==1
            axis equal
        end
    end
end

% Draw the floor
DrawFloor(Sim.Env,Sim.Out.SuppPos(1,1),Sim.Out.SuppPos(end,1));

end

% %%%%%%%% Auxiliary nested functions %%%%%%%% %
% %%%% Draw Circle %%%% %
% Draws a circle of radius R in pos (x,y,z)
function DrawCircle(CB, x, y, z, R, color, ID)
    coordX=zeros(1,CB.CircRes);
    coordY=zeros(1,CB.CircRes);
    coordZ=zeros(1,CB.CircRes);

    for r=1:CB.CircRes
        coordX(1,r) = x+R*cos(r/CB.CircRes*2*pi);
        coordY(1,r) = y+R*sin(r/CB.CircRes*2*pi);
        coordZ(1,r) = z;
    end

    h = patch(coordX,coordY,coordZ,color);
    set(h,'EdgeColor',color.^4);
    set(h,'LineWidth',2*CB.LineWidth);

%         switch ID
%             case 1
%                 CB.RenderObj.Cm1=h;
%             case 2
%                 CB.RenderObj.Cm2=h;
%             case 3
%                 CB.RenderObj.Cmh=h;
%             otherwise
%                 return;
%         end                    
end

% %%%% Draw Link %%%% %
% Draws a link of from (x0,y0) to (x1,y1)
function DrawLink(CB, x0, y0, x1, y1, z)
    Length = sqrt((x1-x0)^2+(y1-y0)^2);
    Center = [(x0+x1)/2;
              (y0+y1)/2];
    Orientation = atan2(y1-y0,x1-x0);

    res.Trans=hgtransform('Parent',gca);

    coordX=zeros(1,2*CB.LinkRes+1);
    coordY=zeros(1,2*CB.LinkRes+1);
    coordZ=zeros(1,2*CB.LinkRes+1);

    x=0;
    y=Length-CB.leg_width/2;

    for r=1:CB.LinkRes
        coordX(1,r) = x+CB.leg_width/2*cos(r/CB.LinkRes*pi);
        coordY(1,r) = y+CB.leg_width/2*sin(r/CB.LinkRes*pi);
        coordZ(1,r) = 0;
    end

    y = CB.leg_width/2;
    for r=CB.LinkRes:2*CB.LinkRes
        coordX(1,r+1) = x+CB.leg_width/2*cos(r/CB.LinkRes*pi);
        coordY(1,r+1) = y+CB.leg_width/2*sin(r/CB.LinkRes*pi);
        coordZ(1,r+1) = 0;
    end

    %         Txy=makehgtform('translate',[Center(1) Center(2) 0]);
    Txy=makehgtform('translate',[x0 y0 z]);
    Rz=makehgtform('zrotate',Orientation-pi/2);
    Sx=makehgtform('scale',[1,Length/CB.L,1]);

    res.Geom=patch(coordX,coordY,coordZ,CB.leg_color);
    set(res.Geom,'EdgeColor',[0 0 0]);
    set(res.Geom,'LineWidth',2*CB.LineWidth);

    set(res.Geom,'Parent',res.Trans);
    set(res.Trans,'Matrix',Txy*Rz*Sx);

end

% %%%% Draw Vector %%%% %
% Draws a vector from x0 to x1
function DrawVector(CB,x0,x1,zIndex,Color)
    VecScale=CB.leg_width*0.75;
    Length=sqrt((x1(1)-x0(1))^2+(x1(2)-x0(2))^2);
    if Length<VecScale*4
        return;
    end

    Dir=(x1-x0)/Length;
    DirPerp=[-Dir(2); Dir(1)];

    Points=zeros(3,7);
    Points(1:2,1)=x0-DirPerp*VecScale/2;
    Points(1:2,2)=x0+DirPerp*VecScale/2;
    Points(1:2,3)=x1-Dir*VecScale*4+DirPerp*VecScale/2;
    Points(1:2,4)=x1-Dir*VecScale*4+DirPerp*VecScale*1.5;
    Points(1:2,5)=x1;
    Points(1:2,6)=x1-Dir*VecScale*4-DirPerp*VecScale*1.5;
    Points(1:2,7)=x1-Dir*VecScale*4-DirPerp*VecScale/2;
    Points(3,:)=zIndex*ones(1,7);

    patch(Points(1,:),Points(2,:),Points(3,:),Color);
end
    
% %%%%%% % Get position % %%%%%% %
function [ x, y ] = GetLegsPos(CB, X, xS, yS, which)
    sS=sin(X(1)); cS=cos(X(1));
    sNS=sin(X(2)); cNS=cos(X(2));

    if strcmp(which,'S')
        x = xS;
        y = yS;
        return;
    end
    if strcmp(which,'m1')
        x = xS-(1-CB.a)*CB.L*sS;
        y = yS+(1-CB.a)*CB.L*cS;
        return;
    end
    if strcmp(which,'Hip')
        x = xS-CB.L*sS;
        y = yS+CB.L*cS;
        return;
    end
    if strcmp(which,'m2')
        x = xS-CB.L*sS+CB.a*CB.L*sNS;
        y = yS+CB.L*cS-CB.a*CB.L*cNS;
        return;
    end
    if strcmp(which,'NS')
        x = xS-CB.L*sS+(CB.L-CB.LegShift)*sNS;
        y = yS+CB.L*cS-(CB.L-CB.LegShift)*cNS;
        return;
    end
    if strcmp(which,'COM')
        % Support mass
        [ m1x, m1y ] = GetLegPos(CB, X, 'm1');
        % Non-support mass
        [ m2x, m2y ] = GetLegPos(CB, X, 'm2');
        % Hip mass
        [ mhx, mhy ] = GetLegPos(CB, X, 'Hip');

        x=(CB.m*m1x+CB.m*m2x+CB.mh*mhx)/(2*CB.m+CB.mh);
        y=(CB.m*m1y+CB.m*m2y+CB.mh*mhy)/(2*CB.m+CB.mh);
    end
end

function DrawFloor(Env,Min,Max)
    FloorX=Min:Env.FloorStep:Max;
    VLStep=(Max-Min)/(Env.VertLines);

    FloorY=Env.Surf(FloorX);

    % Draw horizontal line
    Env.FloorLine=line(FloorX,FloorY, 'LineWidth', 3*Env.LineWidth, 'Color', Env.FloorColor);

    % Draw vertical lines
    for v=1:Env.VertLines
        Env.FloorVLx(v)=Min+v*VLStep;
        Env.FloorVL(v)=line([Env.FloorVLx(v) Env.FloorVLx(v)-1/15],[Env.Surf(Env.FloorVLx(v)) Env.Surf(Env.FloorVLx(v))-1/5],...
                           'LineWidth', 2*Env.LineWidth, 'Color', Env.FloorColor);
    end
end