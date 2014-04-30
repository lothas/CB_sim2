classdef Terrain < handle
    % Version 0.3 - 24/04/2014
    
    % Different terrains: inclined plane, sinusoidal
    %                     infinite parabolla, finite parabolla
    
    properties
        Type=0;
        % 0 - inclined plane
        % 1 - sinusoidal
        % 2 - infinite parabolla
        % 3 - finite parabolla
        
        % Sine terrain params
        sinAmp=0.1;
        sinFreq=1;
        
        % Parabolla constant parK/2*x^2
        parK=0.025 % approx 1.5 deg/m
        
        % Parabolla params
        incline=0;      % Incline (1) up, (-1) down
        start_slope=0;
        end_slope=0;
        
        start_x=0;
        start_y=0;
        end_x=0;
        end_y=0;
        
        % Render parameters
        FloorStep=0.05;
        VertLines=10;
        FloorColor=[0.1,0.3,0];

        FloorLine=0;
        FloorVLx=zeros(1,10);
        FloorVL=zeros(1,10);
        
        LineWidth=1;
        
        % Function handles
        Surfh;
        SurfSlopeh;
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function [Te] = Terrain(varargin)
            switch nargin
                case 0  % create empty terrain (flat)
                    Te; %#ok<VUNUS>
                case 1  % input is terrain type
                    Te.Type=varargin{1};
                case 2  % input is terrain type + const. slope
                    Te.Type=varargin{1};
                    Te.start_slope=varargin{2};
                    Te.end_slope=varargin{2};
                case 3  % input is type + depends on type
                    Te.Type=varargin{1};
                    switch varargin{1}
                        case 1  % Sine amp and fre
                            Te.sinAmp=varargin{2};
                            Te.sinFreq=varargin{3};
                        case 2  % parabolla init slope + curvature
                            Te.start_slope=varargin{2};
                            Te.parK=varargin{3};
                        case 3  % parabolla init and end slope
                            Te.start_slope=varargin{2};
                            Te.end_slope=varargin{3};
                        otherwise
                            Te.start_slope=varargin{2};
                            Te.end_slope=varargin{2};
                    end
                case 4
                    Te.Type=varargin{1};
                    Te.start_slope=varargin{2};
                    Te.end_slope=varargin{3};
                    Te.parK=varargin{4};
            end
            
            % Set init and end points
            Te=SetEndConditions(Te);
            Te.Surfh = {@Te.Surf0, @Te.Surf1, @Te.Surf2, @Te.Surf3};
            Te.SurfSlopeh = {@Te.SurfSlope0, @Te.SurfSlope1,...
                            @Te.SurfSlope2, @Te.SurfSlope3};
        end
        
        function [Te] = SetEndConditions(Te)
            if Te.Type==2
                % For the inifinite parabolla, set the incline
                % based on the curvature sign
                Te.incline=sign(Te.parK);
                Te.parK=Te.parK*Te.incline;
            else
                % Set the inclination based on begin and end slopes
                Te.incline=sign(Te.end_slope-Te.start_slope);
            end
            
            if Te.incline==0
                Te.end_x=0;
                Te.end_y=0;
            else
                % %% This is used for type 3 only (finite par.) %%
                % end_x is actually the distance required
                % from x=start_x to reach the desired angle
                Te.end_x=Te.start_x+(tan(Te.end_slope*pi/180)-tan(Te.start_slope*pi/180))/(Te.incline*Te.parK);
                Te.end_y=Te.start_y+tan(Te.start_slope*pi/180)*(Te.end_x-Te.start_x)+...
                    Te.incline*Te.parK/2*(Te.end_x-Te.start_x)^2;
            end
        end
        
        function [Te] = SetSmoothness(Te,Kslope)
            % The higher the value of Kslope, the faster the slope
            % will change from start_slope to end_slope
            % Default is 1.5 degrees per meter
            Te.parK=Kslope;
            Te=Te.SetEndConditions();
        end
        
        function [y, Trans] = Surf(Te,x)
            if length(x)>1
                ID1 = find(x>=Te.start_x,1,'first');
                ID2 = find(x>=Te.end_x,1,'first');
                if isempty(ID1)
                    ID1 = 1;
                end
                if isempty(ID2)
                    ID2 = length(x);
                end
                [y1, ~] = Te.Surfh{Te.Type+1}(x(1:ID1-1));
                [y2, ~] = Te.Surfh{Te.Type+1}(x(ID1:ID2-1));
                [y3, ~] = Te.Surfh{Te.Type+1}(x(ID2:end));
                y = [y1,y2,y3];
                Trans=[];
            else
                [y,Trans] = Te.Surfh{Te.Type+1}(x);
            end
        end

        function [alpha]=SurfSlope(Te,x)
            alpha = Te.SurfSlopeh{Te.Type+1}(x);
        end
        
% %%%%%%%%%%%% Type 0 - inclined plane %%%%%%%%%%%%
        function [y, Trans] = Surf0(Te,x)
            alpha = Te.SurfSlope(x);
            Trans=[cos(alpha), -sin(alpha);
                   sin(alpha), cos(alpha)];
               
            y=Te.start_y+x*tand(Te.start_slope);
        end

        function [alpha]=SurfSlope0(Te,x) %#ok<INUSD>
            alpha=Te.start_slope*pi/180;
        end
        
% %%%%%%%%%%%% Type 1 - sinusoidal %%%%%%%%%%%%
        function [y, Trans] = Surf1(Te,x)
            alpha = Te.SurfSlope(x);
            Trans=[cos(alpha), -sin(alpha);
                   sin(alpha), cos(alpha)];
               
            if x<Te.start_x
                y=0;            
            else
                y=Te.sinAmp*(1-cos(Te.sinFreq*(x-Te.start_x)));
            end
        end

        function [alpha]=SurfSlope1(Te,x)
            if x<Te.start_x
                alpha=0;
            else
                alpha=Te.sinAmp*Te.sinFreq*sin(Te.sinFreq*(x));
            end
        end
        
% %%%%%%%%%%%% Type 2 - infinite parabolla %%%%%%%%%%%%
        function [y, Trans] = Surf2(Te,x)
            alpha = Te.SurfSlope(x);
            Trans=[cos(alpha), -sin(alpha);
                   sin(alpha), cos(alpha)];
               
            if x<Te.start_x
                y=tan(alpha)*(x-Te.start_x);            
            else
                y=Te.start_y+tan(Te.start_slope*pi/180)*(x-Te.start_x)+Te.incline*Te.parK/2*(x-Te.start_x).^2;
            end
        end

        function [alpha]=SurfSlope2(Te,x)
            if x<Te.start_x
                alpha=Te.start_slope*pi/180;
            else
                alpha=atan(tan(Te.start_slope*pi/180)+Te.incline*Te.parK*(x-Te.start_x));
            end
        end
        
% %%%%%%%%%%%% Type 3 - finite parabolla %%%%%%%%%%%%
        function [y, Trans] = Surf3(Te,x)
            alpha = Te.SurfSlope(x);
            if x<Te.start_x
                y=Te.start_y+tan(Te.start_slope*pi/180)*(x-Te.start_x);            
            else
                if x<Te.end_x
                    y=Te.start_y+tan(Te.start_slope*pi/180)*(x-Te.start_x)+Te.incline*Te.parK/2*(x-Te.start_x).^2;
                else
                    y=Te.end_y+tan(Te.end_slope*pi/180)*(x-Te.end_x);
                end
            end

            Trans=[cos(alpha), -sin(alpha);
                   sin(alpha), cos(alpha)];
        end

        function [alpha]=SurfSlope3(Te,x)
            if x<Te.start_x
                alpha=Te.start_slope*pi/180;
            else
                if x<Te.end_x
                    alpha=atan(tan(Te.start_slope*pi/180)+Te.incline*Te.parK*(x-Te.start_x));
                else
                    alpha=Te.end_slope*pi/180;
                end
            end
        end
        
        % %%%%% Rendering function in Render.m %%%%%
    end
end