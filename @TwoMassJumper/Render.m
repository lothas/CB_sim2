% %%%%%% % Render Two Mass Jumper % %%%%%% %
function TMJ = Render(TMJ,X)
    if isempty(TMJ.RenderObj)
        % Model hasn't been rendered yet
        % Render masses as boxes
        TMJ = DrawBox(TMJ, 0, X(1), 0, 1);
        TMJ = DrawBox(TMJ, 0, X(3), 0, 2);
        
        % Render spring
        bH = TMJ.m_size*1*TMJ.m2;
        TMJ = DrawSpring(TMJ, [0, X(1)], [0, X(3)+bH], 0);
        
        % Render force
%         TMJ = DrawVector(TMJ, [0, X(1)], [0, TMJ.Force], 0, 1);
%         TMJ = DrawVector(TMJ, [0, X(2)], [0, -TMJ.Force], 0, 2);

        % Create transformation handles
        TMJ.RenderObj.Cm1t = hgtransform('Parent',gca);
        TMJ.RenderObj.Cm2t = hgtransform('Parent',gca);
        set(TMJ.RenderObj.Cm1,'Parent',TMJ.RenderObj.Cm1t);
        set(TMJ.RenderObj.Cm2,'Parent',TMJ.RenderObj.Cm2t);

        % Save initial position as reference
        TMJ.OldmPos=[X(1), X(3)];

        % Render additional parameters (m, L, I)
%         text(x1+0.02*TMJ.L, y1+0.2*TMJ.L, '{\itm_h}','FontSize',24,'FontWeight','bold','HorizontalAlignment','Center');
%         text(x1a-0.04*TMJ.L, y1a+0.12*TMJ.L, '{\itm, L, I}','FontSize',24,'FontWeight','bold','HorizontalAlignment','Right');

        axis equal
        
        % Finished rendering
        % Call function again to proceed with the code below
        TMJ = TMJ.Render(X);
    else
        if ishandle(TMJ.RenderObj.Cm1)==0
            % If window was closed and handle destroyed, re-render
            TMJ.RenderObj=[];
            TMJ = TMJ.Render(X);
        else
            % Model was already rendered
                        
            % Use transformations to translate the masses
            Tm1xy = makehgtform('translate',[0 X(1)-TMJ.OldmPos(1) 0]);
            Tm2xy = makehgtform('translate',[0 X(3)-TMJ.OldmPos(2) 0]);
            set(TMJ.RenderObj.Cm1t,'Matrix',Tm1xy);
            set(TMJ.RenderObj.Cm2t,'Matrix',Tm2xy);
        end
    end

    % %%%%%%%% Auxiliary nested functions %%%%%%%% %
    % %%%% Draw Box %%%% %
    % Draws a box of size m*(2x1) in pos (x,y,z)
    function [ TMJ ] = DrawBox(TMJ, x, y, z, ID)
        switch ID
            case 1
                W_2 = TMJ.m_size*2*TMJ.m1/2;
                H = TMJ.m_size*1*TMJ.m1;
            case 2
                W_2 = TMJ.m_size*2*TMJ.m2/2;
                H = TMJ.m_size*1*TMJ.m2;
            otherwise
                return;
        end            
        
        coordX = [x - W_2; x - W_2; x + W_2; x + W_2];
        coordY = [y; y + H; y + H; y];
        coordZ = [z; z; z; z];
        
        h = patch(coordX,coordY,coordZ,TMJ.m_color);
        set(h,'EdgeColor',TMJ.m_color.^4);
        set(h,'LineWidth',TMJ.LineWidth);

        switch ID
            case 1
                TMJ.RenderObj.Cm1 = h;
            case 2
                TMJ.RenderObj.Cm2 = h;
        end       
    end

    % %%%% Draw Spring %%%% %
    % Draws a spring from x0 to x1
    function [ TMJ ] = DrawSpring(TMJ, PosS, PosE, Z)
        NumTurns=5;

        Center = (PosS+PosE)/2;
        Length = sqrt((PosE(1)-PosS(1))^2+(PosE(2)-PosS(2))^2);
        Orientation = atan2(PosE(2)-PosS(2),PosE(1)-PosS(1));

        coordX = zeros(1,NumTurns*2+2);
        coordY = zeros(1,NumTurns*2+2);
        Z = ones(1,NumTurns*2+2)*Z;

        Step = Length/NumTurns;

        coordX(1) = -Length/2;
        coordX(end) = Length/2;
        for t=1:NumTurns
            coordX(2*t)=coordX(1)+Step*(1/4+t-1);
            coordX(2*t+1)=coordX(1)+Step*(3/4+t-1);
            coordY(2*t)=TMJ.m_size;
            coordY(2*t+1)=-TMJ.m_size;
        end

        % Rotate spring
        Coords = [coordX; coordY];
        RotMatrix = [cos(Orientation) -sin(Orientation);
                   sin(Orientation)  cos(Orientation)];
        Coords = repmat(Center,1,NumTurns*2+2)+RotMatrix*Coords;

        if isempty(TMJ.RenderObj)
            TMJ.RenderObj.Spring = ...
                plot3(Coords(1,:), Coords(2,:), Z, 'Color', [0, 0, 0], 'LineWidth', 2*TMJ.LineWidth);
        else
            set(TMJ.RenderObj.Spring(1),'XData',Coords(1,:));
            set(TMJ.RenderObj.Spring(1),'YData',Coords(2,:));
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