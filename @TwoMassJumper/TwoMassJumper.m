classdef TwoMassJumper
    % Version 0.1 - 02/06/2014
    
    % Linear model of two masses connected by a spring.
    % Upper mass m1, lower mass m2.
    
    properties
        % %%%%%% % System parameters % %%%%%% %
        m1 = 2; % top mass
        m2 = 0.2; % bottom mass
        
        k = 1; % spring constant
        l0 = 0.5; % spring zero length
        
        g = 9.81; % Gravity
        
        % Dampening
        damp = 0;
        
        % enum
        Upper = 1;
        Lower = 2;
        
        % Nondimensional Parameters
        M = 0; omega = 0; dist = 0;
        kappa_c = 0; gamma = 0;
        eta=0;
        
        % ND matrices elements
        MM1; M1;
        
        % External force
        Force = 0;
        
        % Runtime variables
        Phase = 'stance'; % or 'flight'
        last_t = 0; % Last time since impact
        
        stDim = 4; % state dimension
        nEvents = 2; % num. of simulation events
        % 1 - Touch down
        % 2 - Lift-off
        
        % Set keys
        SetKeys = {'m1','m2','k','l0','g','damp','Phase',...
            'm_size','m_color','LineWidth'};
                
        % %%%%%% % Render parameters % %%%%%% %
        m_size = 0.1;
        m_color=[0.2,0.6,0.8];

        OldmPos=0;
        
        RenderObj;
        RenderVectors=0;
        RenderStance=0;
        RenderParams=0;
        
        LineWidth=2;
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function TMJ = TwoMassJumper(varargin)
            switch nargin
                case 0 % Empty object
                    TMJ; %#ok<VUNUS>
                case 4 % Input is: m1, m2, k, damp
                    TMJ.m1 = varargin{1}; % top mass
                    TMJ.m2 = varargin{2}; % bottom mass
                    TMJ.k = varargin{3}; % spring constant
                    TMJ.damp = varargin{4}; % damping

                    TMJ.g=9.81;
                case 5 % Input is: m1, m2, k, l0, damp
                    TMJ.m1 = varargin{1}; % leg mass
                    TMJ.m2 = varargin{2}; % hip mass
                    TMJ.k = varargin{3}; % leg length
                    TMJ.l0 = varargin{4}; % leg center of mass
                    TMJ.damp = varargin{5}; % damping

                    TMJ.g=9.81;                
            end
            TMJ = TMJ.SetND();
        end
        
        function [TMJ] = SetND(TMJ)
            % Nondimensional Parameters
            TMJ.M = TMJ.m2/TMJ.m1;
            TMJ.omega = sqrt(TMJ.k*(TMJ.m2+TMJ.m1)/TMJ.m2/TMJ.m1);
            TMJ.dist = TMJ.g*TMJ.m2/TMJ.k;
            % Coeff that multiplies the force:
            TMJ.kappa_c = 1/TMJ.g/TMJ.m2;
            TMJ.gamma = TMJ.l0*TMJ.k/TMJ.g/TMJ.m2;
            TMJ.eta = TMJ.damp*sqrt(1/TMJ.m1+1/TMJ.m2)/sqrt(TMJ.k);
            
            % ND matrices elements
            TMJ.M1 = 1/(TMJ.M+1);
            TMJ.MM1 = TMJ.M*TMJ.M1;
        end
        
        % %%%%%% % Get position % %%%%%% %
        function [ y ] = GetPos(TMJ, X, which)
            switch which
                case 'm1'
                    y = X(:,1);
                case 'm2'
                    y = X(:,3);
                case 'COM'
                    % Upper mass
                    m1y = TMJ.GetPos(X, 'm1');
                    % Lower mass
                    m2y = TMJ.GetPos(X, 'm2');

                    y = (TMJ.m1*m1y+TMJ.m2*m2y)/(TMJ.m1+TMJ.m2);
            end
        end
        
        % %%%%%% % Get velocity % %%%%%% %
        function [ ydot ] = GetVel(TMJ, X, which)
            switch which
                case 'm1'
                    ydot = X(:,2);
                case 'm2'
                    ydot = X(:,4);
                case 'COM'
                    % Upper mass
                    m1y = TMJ.GetVel(X, 'm1');
                    % Lower mass
                    m2y = TMJ.GetVel(X, 'm2');

                    ydot = (TMJ.m1*m1y+TMJ.m2*m2y)/(TMJ.m1+TMJ.m2);
            end
        end
        
        function [XND] = D2ND(TMJ,X)
            % Transform coordiantes to non-dimensional
            XND = [X(:,1)/TMJ.dist,...
                   X(:,2)/TMJ.dist/TMJ.omega,...
                   X(:,3)/TMJ.dist,...
                   X(:,4)/TMJ.dist/TMJ.omega];
        end
        
        function [X] = ND2D(TMJ,XND)
            % Transform non-dimensional coordiantes to regular coord.
            X = [XND(:,1)*TMJ.dist,...
                 XND(:,2)*TMJ.dist*TMJ.omega,...
                 XND(:,3)*TMJ.dist,...
                 XND(:,4)*TMJ.dist*TMJ.omega];
        end
        
        % %%%%%% % Get ground reaction forces % %%%%%% %
        function [F] = GetGRF(TMJ, X)
            % Calculate GRF based on the 2 DOF model (y1,y2)
            [XND] = TMJ.D2ND(X);
            
            F = TMJ.gama + TMJ.kappa_c*TMJ.Force - XND(:,1) - TMJ.eta*XND(:,2) + 1;
        end
        
        function [Weight] = GetWeight(TMJ)
            Weight = (TMJ.m1+TMJ.m2)*TMJ.g;
        end
        
        function [KE] = GetKineticEnergy(TMJ, X)
            m1y = TMJ.GetVel(X, 'm1');
            m2y = TMJ.GetVel(X, 'm2');
            KE = 1/2*TMJ.m1*m1y.^2+1/2*TMJ.m2*m2y.^2;
        end
        
        function [PE] = GetPotentialEnergy(TMJ, X)
            PE = TMJ.g*TMJ.m1*X(:,1) + TMJ.g*TMJ.m2*X(:,3) + ...
                 1/2*TMJ.k*(X(:,1)-X(:,2)-TMJ.l0).^2;
        end
        
        function [A, Bu] = GetSystemMat(TMJ)
            switch TMJ.Phase
                case 'stance'
                    kappa = TMJ.kappa_c*TMJ.Force;
                    A = [0 1 0 0;
                        -TMJ.MM1 -TMJ.MM1*TMJ.eta 0 0;
                         0 0 0 0;
                         0 0 0 0];
                    Bu = [0; TMJ.MM1*(TMJ.gamma+kappa)-TMJ.M1; 0; 0];
                case 'flight'
                    A = [0 1 0 0;
                        -TMJ.MM1 -TMJ.MM1*TMJ.eta TMJ.MM1 TMJ.MM1*TMJ.eta;
                         0 0 0 1;
                         TMJ.M1 TMJ.M1*TMJ.eta -TMJ.M1 -TMJ.M1*TMJ.eta];
                    Bu = [0;
                          TMJ.MM1*TMJ.gamma-TMJ.M1;
                          0;
                          -TMJ.M1*(TMJ.gamma+1)];
            end
        end
            
        % %%%%%% % Derivative % %%%%%% %
        function [Xdot] = Derivative(TMJ, t, X) %#ok<INUSL>
            % Derivative uses normal coordinates
            % Translate first to nondimensional
            XND = TMJ.D2ND(X);

            [A, Bu] = TMJ.GetSystemMat();
            XNDdot = A*XND'+Bu;
            
            % Return to normal coordinates
            Xdot = TMJ.ND2D(XNDdot);
        end
        
        % %%%%%% % Calculate Impact % %%%%%% %
        function [Xa] = CalcImpact(TMJ, Xb) %#ok<MANU>
            % Calculate Impact
            Xa = [Xb(1); Xb(2); Xb(3); 0];
            
            if abs(Xb(3)>1e-4)
                disp('ERROR: Impact occurred when height isn''t 0');
                disp(['Height = ',num2str(Xb(3))]);
            end
        end
        
        % %%%%%% % Events % %%%%%% %
        function [value, isterminal, direction] = Events(TMJ, X, Floor)
            value=ones(TMJ.nEvents,1);
            isterminal=ones(TMJ.nEvents,1);
            direction=-ones(TMJ.nEvents,1);
            
            switch TMJ.Phase
                case 'stance'
                    value(2) = -TMJ.GetGRF(X);
                case 'flight'
                    y = TMJ.GetPos(X,'m2');
                    value(1) = y - Floor.Surf(0);
            end
        end
        
        % %%%%%% % Events % %%%%%% %
        function [TMJ,Xa] = HandleEvent(TMJ, evID, Xb, t) %#ok<INUSD>
            switch evID
                case 1
                    % Touch down
                    % Calculate Impact
                    Xa = TMJ.CalcImpact(Xb);

                    % Set phase to stance
                    TMJ.Phase = 'stance';
                case 2
                    % Lift off
                    Xa = Xb;
                    % Set phase to flight
                    TMJ.Phase = 'flight';
            end
        end
                
        % %%%%%% % Render Compass Biped % %%%%%% %
        % In Render.m
        % Includes also DrawCircle, DrawLink and DrawVector
        
    end % end methods
end