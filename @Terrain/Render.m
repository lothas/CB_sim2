% %%%%%%%%%%%% Rendering function %%%%%%%%%%%%
function [Te]=Render(Te,Min,Max)
    FloorX=Min:Te.FloorStep:Max;
    VLStep=(Max-Min)/(Te.VertLines);

    FloorY=Te.Surf(FloorX);

    if Te.FloorLine==0 || ishandle(Te.FloorLine)==0
        % Draw horizontal line
        Te.FloorLine=line(FloorX,FloorY, 'LineWidth', 3*Te.LineWidth, 'Color', Te.FloorColor);

        % Draw vertical lines
        for v=1:Te.VertLines
            Te.FloorVLx(v)=Min+v*VLStep;
            Te.FloorVL(v)=line([Te.FloorVLx(v) Te.FloorVLx(v)-1/15],[Te.Surf(Te.FloorVLx(v)) Te.Surf(Te.FloorVLx(v))-1/5],...
                               'LineWidth', 2*Te.LineWidth, 'Color', Te.FloorColor);
        end
    else
        % Update horizontal line
        set(Te.FloorLine, 'XData', FloorX);
        set(Te.FloorLine, 'YData', FloorY);

        % Update vertical lines
        if Te.FloorVLx(1)<Min
            Te.FloorVLx(1:end-1)=Te.FloorVLx(2:end);
            Te.FloorVLx(end)=Te.FloorVLx(end-1)+VLStep;

            for v=1:Te.VertLines
                set(Te.FloorVL(v), 'XData', [Te.FloorVLx(v) Te.FloorVLx(v)-1/15]);
                set(Te.FloorVL(v), 'YData', [Te.Surf(Te.FloorVLx(v)) Te.Surf(Te.FloorVLx(v))-1/5]);
            end
        end

        if Te.FloorVLx(end)>Max
            Te.FloorVLx(2:end)=Te.FloorVLx(1:end-1);
            Te.FloorVLx(1)=Te.FloorVLx(2)-VLStep;

            for v=1:Te.VertLines
                set(Te.FloorVL(v), 'XData', [Te.FloorVLx(v) Te.FloorVLx(v)-1/15]);
                set(Te.FloorVL(v), 'YData', [Te.Surf(Te.FloorVLx(v)) Te.Surf(Te.FloorVLx(v))-1/5]);
            end
        end
    end        
end