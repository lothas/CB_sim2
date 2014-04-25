% %%%%%% % Render Compass Biped % %%%%%% %
function CB = Render(CB,X)
    % Get the position of each mass and the non-support leg
    % Support leg position is (CB.xS, CB.yS)
    [x1,y1]=CB.GetPos(X,'Hip');
    [x2,y2]=CB.GetPos(X,'NS');
    [x1a,y1a]=CB.GetPos(X,'m1');
    [x2a,y2a]=CB.GetPos(X,'m2');

    if isempty(CB.RenderObj)
        % Model hasn't been rendered yet
        % Render masses as circles
        CB=DrawCircle(CB, x1a, y1a, 0, CB.m_radius, CB.m_color,1);
        CB=DrawCircle(CB, x2a, y2a, 0, CB.m_radius, CB.m_color,2);
        CB=DrawCircle(CB, x1, y1, 0, CB.mh_radius, CB.mh_color,3);
        % Render links
        CB.RenderObj.nL1=DrawLink(CB, CB.xS, CB.yS, x1, y1, 0, []);
        CB.RenderObj.nL2=DrawLink(CB, x1, y1, x2, y2, 0, []);

        % Create transformation handles
        CB.RenderObj.Cm1t = hgtransform('Parent',gca);
        CB.RenderObj.Cm2t = hgtransform('Parent',gca);
        CB.RenderObj.Cmht = hgtransform('Parent',gca);
        set(CB.RenderObj.Cm1,'Parent',CB.RenderObj.Cm1t);
        set(CB.RenderObj.Cm2,'Parent',CB.RenderObj.Cm2t);
        set(CB.RenderObj.Cmh,'Parent',CB.RenderObj.Cmht);

        % Save initial position as reference
        CB.OldmPos=[x1a, y1a, x2a, y2a, x1, y1];

        % Draw vectors as required
        % (These are used for calculating the
        % conservation of angular momentum)
        switch CB.RenderVectors
            case 1 % Whole body around new stance foot before
                % Draw 'r' (lever) vectors
                DrawVector(CB, [x2;y2],[x1;y1],5,[1 0 0]);
                DrawVector(CB, [x2;y2],[x1a;y1a],5,[1 0 0]);
                DrawVector(CB, [x2;y2],[x2a;y2a],6,[1 0 0]);

                % Draw 'V' (velocity) vectors
                [Vxmh,Vymh]=CB.GetVel(X,'Hip');
                [Vxm1,Vym1]=CB.GetVel(X,'m1');
                [Vxm2,Vym2]=CB.GetVel(X,'m2');
                DrawVector(CB, [x1a;y1a],[x1a+Vxm1;y1a+Vym1],5,[0 0.5 1]);
                DrawVector(CB, [x2a;y2a],[x2a+Vxm2;y2a+Vym2],5,[0 0.5 1]);
                DrawVector(CB, [x1;y1],[x1+Vxmh;y1+Vymh],5,[0 0.5 1]);
            case 2 % Whole body around new stance foot after
                % Draw 'r' (lever) vectors
                DrawVector(CB, [CB.xS;CB.yS],[x1;y1],5,[1 0 0]);
                DrawVector(CB, [CB.xS;CB.yS],[x1a;y1a],6,[1 0 0]);
                DrawVector(CB, [CB.xS;CB.yS],[x2a;y2a],5,[1 0 0]);

                % Draw 'V' (velocity) vectors
                [Vxmh,Vymh]=CB.GetVel(X,'Hip');
                [Vxm1,Vym1]=CB.GetVel(X,'m1');
                [Vxm2,Vym2]=CB.GetVel(X,'m2');
                DrawVector(CB, [x1a;y1a],[x1a+Vxm1;y1a+Vym1],5,[0 0.5 1]);
                DrawVector(CB, [x2a;y2a],[x2a+Vxm2;y2a+Vym2],5,[0 0.5 1]);
                DrawVector(CB, [x1;y1],[x1+Vxmh;y1+Vymh],5,[0 0.5 1]);
            case 3 % New swing leg around hip before
                % Draw 'r' (lever) vectors
                DrawVector(CB, [x1;y1],[x1a;y1a],5,[1 0 0]);

                % Draw 'V' (velocity) vectors
                [Vxm1,Vym1]=CB.GetVel(X,'m1');
                DrawVector(CB, [x1a;y1a],[x1a+Vxm1;y1a+Vym1],5,[0 0.5 1]);
            case 4 % New swing leg around hip after
                % Draw 'r' (lever) vectors
                DrawVector(CB, [x1;y1],[x2a;y2a],5,[1 0 0]);

                % Draw 'V' (velocity) vectors
                [Vxm2,Vym2]=CB.GetVel(X,'m2');
                DrawVector(CB, [x2a;y2a],[x2a+Vxm2;y2a+Vym2],5,[0 0.5 1]);
            otherwise
        end

        % Add some lines to show foot impact/detachment from ground
        switch CB.RenderStance
            case 1 % Leg impacts ground
                % Show "hinge" around stance foot
                DrawCircle(CB, CB.xS, CB.yS, 5, CB.leg_width/2, [0,0,0], 0);

                % Draw lines around swing foot
                pLineWidth=3*CB.LineWidth;
                LineColor=[0.8, 0.5, 0];
                Aperture=90;
                for i=1:5
                    Angle=(i-3)/4*Aperture-(X(2)-X(1))*180/pi;
                    PointsX=[x2-CB.leg_width*sind(Angle) x2-5*CB.leg_width*sind(Angle)];
                    PointsY=[y2-CB.leg_width*cosd(Angle) y2-5*CB.leg_width*cosd(Angle)];
                    line(PointsX,PointsY,'LineWidth',pLineWidth,'Color',LineColor);
                end
            case 2 % Leg detaches from ground
                % Gotta change this part, it's the same as above...
                DrawCircle(CB, CB.xS, CB.yS, 5, CB.leg_width/2, [0,0,0], 0);

                pLineWidth=3*CB.LineWidth;
                LineColor=[0.8, 0.5, 0];
                Aperture=90;
                for i=1:5
                    Angle=(i-3)/4*Aperture+(X(2)-X(1))*180/pi;
                    PointsX=[x2-CB.leg_width*sind(Angle) x2-5*CB.leg_width*sind(Angle)];
                    PointsY=[y2+CB.leg_width*cosd(Angle) y2+5*CB.leg_width*cosd(Angle)];
                    line(PointsX,PointsY,'LineWidth',pLineWidth,'Color',LineColor);
                end
            otherwise                        
        end

        % Render params (m, L and angles)
        if CB.RenderParams
            MinAngle=10*pi/180;
            
            % Draw arc for stance leg angle, theta 1
            if abs(X(1))>MinAngle
                % Draw arc
                NewL=(0.75-abs(X(1))/(2*pi)*0.5)*CB.L;
                N=ceil(abs(X(1))/(2*pi)*180*(NewL/CB.L/0.75));
                Points1=zeros(2,N);
                TempL=CB.L;
                CB.L=NewL;
                Xpath=X;
                for i=1:N
                    Xpath(1)=(i-1)/N*X(1);
                    [x, y]=CB.GetPos(Xpath,'Hip');
                    Points1(:,i)=[x; y];
                    % Points1(:,i)=[CB.xS+0.75*CB.L*sin((i-1)/N*X(1)); CB.yS+0.75*CB.L*cos((i-1)/N*X(1))];
                end   
                for i=1:floor(N/2)-1
                    line([Points1(1,2*i-1),Points1(1,2*i)],[Points1(2,2*i-1),Points1(2,2*i)],'LineStyle','--','LineWidth',3*CB.LineWidth,'Color',[0 0 0]);
                end
                CB.L=TempL;

                % Draw vertical line
                line([CB.xS CB.xS],[CB.yS CB.yS+NewL+0.1*CB.L],'LineStyle','--','LineWidth',3*CB.LineWidth,'Color',[0 0 0]);

                % Draw arrow
                VecScale=CB.leg_width;
                v0=[Points1(1,end-3); Points1(2,end-3)];
                v1=[Points1(1,end); Points1(2,end)];
                Length=sqrt((v1(1)-v0(1))^2+(v1(2)-v0(2))^2);

                Dir=(v1-v0)/Length;
                DirPerp=[-Dir(2); Dir(1)];

                arPoints1=zeros(3,3);
                arPoints1(1:2,1)=v1-Dir*3*VecScale-DirPerp*VecScale;
                arPoints1(1:2,2)=v1-Dir*3*VecScale+DirPerp*VecScale;
                arPoints1(1:2,3)=v1;
                arPoints1(3,:)=5*ones(1,3);

                patch(arPoints1(1,:),arPoints1(2,:),arPoints1(3,:),[0 0 0]);

                % Print label
                if abs(X(1))<2*MinAngle
                    if X(1)>0
                        text(Points1(1,1)+3*CB.leg_width,Points1(2,1),'\theta_1','FontSize',24,'FontWeight','bold');
                    else
                        text(Points1(1,1)-3*CB.leg_width,Points1(2,1),'\theta_1','FontSize',24,'FontWeight','bold');
                    end
                else
                    % Print label outside angle arc
                    TempL=CB.L;
                    CB.L=NewL+0.04*TempL;
                    Xpath=X;
                    if abs(X(1))>pi
                        Xpath(1)=0.9*X(1);
                    else
                        Xpath(1)=X(1)/2;
                    end
                    [x, y]=CB.GetPos(Xpath,'Hip');
                    CB.L=TempL;
                    text(x,y,'\theta_1','FontSize',24,'FontWeight','bold');                            
                end
            end
            
            % Draw arc for swing leg angle, theta 2
            if abs(X(2))>MinAngle
                [BaseX, BaseY]=CB.GetPos(X,'NS');

                % Draw arc
                NewL=(0.75-abs(X(2))/(2*pi)*0.5)*CB.L;
                N=ceil(abs(X(2))/(2*pi)*180*(NewL/CB.L/0.75));
                Points2=zeros(2,N);
                for i=1:N
                    NewTheta=(i-1)/N*X(2);
                    Points2(:,i)=[BaseX-NewL*sin(NewTheta); BaseY+NewL*cos(NewTheta)];
                end   
                for i=1:floor(N/2)-1
                    line([Points2(1,2*i-1),Points2(1,2*i)],[Points2(2,2*i-1),Points2(2,2*i)],'LineStyle','--','LineWidth',3*CB.LineWidth,'Color',[0 0 0]);
                end

                % Draw vertical line
                line([BaseX BaseX],[BaseY BaseY+NewL+0.1*CB.L],'LineStyle','--','LineWidth',3*CB.LineWidth,'Color',[0 0 0]);

                % Draw arrow
                VecScale=CB.leg_width;
                v0=[Points2(1,end-3); Points2(2,end-3)];
                v1=[Points2(1,end); Points2(2,end)];
                Length=sqrt((v1(1)-v0(1))^2+(v1(2)-v0(2))^2);

                Dir=(v1-v0)/Length;
                DirPerp=[-Dir(2); Dir(1)];

                arPoints2=zeros(3,3);
                arPoints2(1:2,1)=v1-Dir*3*VecScale-DirPerp*VecScale;
                arPoints2(1:2,2)=v1-Dir*3*VecScale+DirPerp*VecScale;
                arPoints2(1:2,3)=v1;
                arPoints2(3,:)=5*ones(1,3);

                patch(arPoints2(1,:),arPoints2(2,:),arPoints2(3,:),[0 0 0]);

                % Print label
                if abs(X(2))<2*MinAngle
                    if X(2)>0
                        text(Points2(1,1)+3*CB.leg_width,Points2(2,1),'\theta_2','FontSize',24,'FontWeight','bold');
                    else
                        text(Points2(1,1)-3*CB.leg_width,Points2(2,1),'\theta_2','FontSize',24,'FontWeight','bold');
                    end
                else
                    % Print label outside angle arc
                    NewL=NewL+0.12*TempL;
                    if abs(X(2))>pi
                        NewTheta=X(2)/7;
                    else
                        NewTheta=X(2)/2;
                    end
                    text(BaseX-NewL*sin(NewTheta), BaseY+NewL*cos(NewTheta),'\theta_2','FontSize',24,'FontWeight','bold');                            
                end
            end

            % Render additional parameters (m, L, I)
            text(x1+0.02*CB.L, y1+0.2*CB.L, '{\itm_h}','FontSize',24,'FontWeight','bold','HorizontalAlignment','Center');
            text(x1a-0.04*CB.L, y1a+0.12*CB.L, '{\itm, L, I}','FontSize',24,'FontWeight','bold','HorizontalAlignment','Right');
        end

        axis equal
        
        % Finished rendering
        % Call function again to proceed with the code below
        Render(CB,X);
    else
        if ishandle(CB.RenderObj.Cm1)==0
            % If window was closed and handle destroyed, re-render
            CB.RenderObj=[];
            Render(CB,X);
        else
            % Model was already rendered
            % Update the objects
            if CB.Support==CB.Left
                Z = -1;
            else
                Z = 1;
            end
            
            % Re-draw links
            DrawLink(CB, CB.xS, CB.yS, x1, y1, Z, CB.RenderObj.nL1);
            DrawLink(CB, x1, y1, x2, y2, -Z, CB.RenderObj.nL2);
            
            % Use transformations to translate the masses
            Tm1xy = makehgtform('translate',[x1a-CB.OldmPos(1) y1a-CB.OldmPos(2) 0.5+Z]);
            Tm2xy = makehgtform('translate',[x2a-CB.OldmPos(3) y2a-CB.OldmPos(4) 0.5-Z]);
            Tmhxy = makehgtform('translate',[x1-CB.OldmPos(5) y1-CB.OldmPos(6) 0]);
            set(CB.RenderObj.Cm1t,'Matrix',Tm1xy);
            set(CB.RenderObj.Cm2t,'Matrix',Tm2xy);
            set(CB.RenderObj.Cmht,'Matrix',Tmhxy);
        end
    end

    % %%%%%%%% Auxiliary nested functions %%%%%%%% %
    % %%%% Draw Circle %%%% %
    % Draws a circle of radius R in pos (x,y,z)
    function [ CB ] = DrawCircle(CB, x, y, z, R, color, ID)
        coordX=zeros(1,CB.CircRes);
        coordY=zeros(1,CB.CircRes);
        coordZ=zeros(1,CB.CircRes);

        for r=1:CB.CircRes
            coordX(1,r)=x+R*cos(r/CB.CircRes*2*pi);
            coordY(1,r)=y+R*sin(r/CB.CircRes*2*pi);
            coordZ(1,r)=z;
        end

        h=patch(coordX,coordY,coordZ,color);
        set(h,'EdgeColor',color.^4);
        set(h,'LineWidth',2*CB.LineWidth);

        switch ID
            case 1
                CB.RenderObj.Cm1=h;
            case 2
                CB.RenderObj.Cm2=h;
            case 3
                CB.RenderObj.Cmh=h;
            otherwise
                return;
        end                    
    end

    % %%%% Draw Link %%%% %
    % Draws a link of from (x0,y0) to (x1,y1)
    function [ res ] = DrawLink(CB, x0, y0, x1, y1, z, Obj)
        if isempty(Obj)
            Length=sqrt((x1-x0)^2+(y1-y0)^2);
            Center=[(x0+x1)/2;
                    (y0+y1)/2];
            Orientation=atan2(y1-y0,x1-x0);

            res.Trans=hgtransform('Parent',gca);
            Txy=makehgtform('translate',[Center(1) Center(2) 0]);
            Rz=makehgtform('zrotate',Orientation-pi/2);

            coordX=zeros(1,2*CB.LinkRes+1);
            coordY=zeros(1,2*CB.LinkRes+1);
            coordZ=zeros(1,2*CB.LinkRes+1);

            x=0;
            y=Length/2-CB.leg_width/2;
            for r=1:CB.LinkRes
                coordX(1,r)=x+CB.leg_width/2*cos(r/CB.LinkRes*pi);
                coordY(1,r)=y+CB.leg_width/2*sin(r/CB.LinkRes*pi);
                coordZ(1,r)=0;
            end

            y=-Length/2+CB.leg_width/2;
            for r=CB.LinkRes:2*CB.LinkRes
                coordX(1,r+1)=x+CB.leg_width/2*cos(r/CB.LinkRes*pi);
                coordY(1,r+1)=y+CB.leg_width/2*sin(r/CB.LinkRes*pi);
                coordZ(1,r+1)=0;
            end

            res.Geom=patch(coordX,coordY,coordZ,CB.leg_color);
            set(res.Geom,'EdgeColor',[0 0 0]);
            set(res.Geom,'LineWidth',2*CB.LineWidth);

            set(res.Geom,'Parent',res.Trans);
            set(res.Trans,'Matrix',Txy*Rz);
        else
            Center=[(x0+x1)/2;
                    (y0+y1)/2];
            Orientation=atan2(y1-y0,x1-x0);
            Length=sqrt((x1-x0)^2+(y1-y0)^2);

            Txy=makehgtform('translate',[Center(1) Center(2) z]);
            Rz=makehgtform('zrotate',Orientation-pi/2);
            Sx=makehgtform('scale',[Length/CB.L,1,1]);
            set(Obj.Trans,'Matrix',Txy*Rz*Sx);
            res=1;
        end
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
end