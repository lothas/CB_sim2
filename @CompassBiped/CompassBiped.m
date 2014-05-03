classdef CompassBiped
    % Version 0.5 - 25/04/2014
    
    % A compass biped. Two kneeless legs connected at the hip.
    % Each leg has a mass "m" at a distance "a" from the hip.
    % Mass can be evenly distributed by setting I=1/3*m*L^2.
    % Hip mass = mh
    
    % Model with theta1 = ccw angle of stance leg from vertical
    %            theta2 = ccw angle of swing leg from vertical
    
    properties
        % %%%%%% % System parameters % %%%%%% %
        m=3; % leg mass
        mh=10; % hip mass
        
        L=1; % leg length
        a=0.5; % leg center of mass
        I=1/3*3*1^2; % leg moment of inertia
        
        LegShift=0; % leg shift var. for clearance
        Clearance=0.03; % clearance

        g=9.81; % Gravity
        
        % Dampening
        dampNS=0; %0.3;
        dampS=0; %0.3;
        
        % enum
        Left=1;
        Right=2;
        
        % Nondimensional Parameters
        alpha=0; aalpha=0;
        omega=0; Igag=0;
        C1gag=0; C2gag=0;
        Ek=0;
        
        % Matrix elements
        M1=0; M2=0; P1=0; P2=0;
        
        % Support Leg position
        xS=0; yS=0;
        
        % Support
        Support=2; % Right
        
        % Control torques
        Torques=[0;0]; % Ankle; Hip
        
        % Runtime variables
        last_t = 0; % Last time since impact
        curSpeed = 0; % Current walking speed
        
        stDim=4; % state dimension
        nEvents=2; % num. of simulation events
        % 1 - Foot contact
        % 2 - Robot falling
        
        % %%%%%% % Render parameters % %%%%%% %
        m_radius=0.015*3;
        m_color=[0.2,0.6,0.8];

        mh_radius=0.015*6;
        mh_color=[0.5,0.5,0.7];
        
        leg_width=0.025;
        leg_color=[0.1, 0.3, 0.8];

        OldmPos=0;
        
        RenderObj;
        RenderVectors=0;
        RenderStance=0;
        RenderParams=0;
        
        CircRes=12;
        LinkRes=10;
        LineWidth=1;
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function CB=CompassBiped(varargin)
            switch nargin
                case 0 % Empty object
                    CB; %#ok<VUNUS>
                case 4 % Input is: m, mh, L, a
                    CB.m=varargin{1}; % leg mass
                    CB.mh=varargin{2}; % hip mass
                    CB.L=varargin{3}; % leg length
                    CB.a=varargin{4}; % leg center of mass

                    CB.I=varargin{1}*varargin{4}^2*varargin{3}^2; % moment of inertia

                    CB.g=9.81;

                    % Support Leg position
                    CB.xS=0;
                    CB.yS=0;

                    % %%%%%% % Render parameters % %%%%%% %
                    CB.m_radius=0.03*varargin{1};
                    CB.mh_radius=0.03*varargin{2};
                case 6 % Input is: m, mh, L, a, x0, y0
                    CB.m=varargin{1}; % leg mass
                    CB.mh=varargin{2}; % hip mass
                    CB.L=varargin{3}; % leg length
                    CB.a=varargin{4}; % leg center of mass

                    CB.I=varargin{1}*varargin{4}^2*varargin{3}^2; % moment of inertia

                    CB.g=9.81;

                    % Support Leg position
                    CB.xS=varargin{5};
                    CB.yS=varargin{6};

                    % %%%%%% % Render parameters % %%%%%% %
                    CB.m_radius=0.03*varargin{1};
                    CB.mh_radius=0.03*varargin{2};
                case 8 % Input is: m, mh, L, a, I, g, x0, y0
                    CB.m=varargin{1}; % leg mass
                    CB.mh=varargin{2}; % hip mass
                    CB.L=varargin{3}; % leg length
                    CB.a=varargin{4}; % leg center of mass
                    CB.I=varargin{5}; % moment of inertia
                    CB.g=varargin{6}; % gravity

                    % Support Leg position
                    CB.xS=varargin{7};
                    CB.yS=varargin{8};

                    % %%%%%% % Render parameters % %%%%%% %
                    CB.m_radius=0.03*varargin{1};
                    CB.mh_radius=0.03*varargin{2};                    
            end
            CB=CB.SetND();
        end
        
        function [CB] = SetND(CB)
            % Nondimensional Parameters
            CB.alpha=CB.m/CB.mh;
            CB.aalpha=CB.a*CB.alpha;
            CB.omega=sqrt(CB.g/CB.L);
            CB.Igag=CB.I/CB.mh/CB.L^2;
            LMg = CB.mh*CB.L*CB.g;
            DampND = CB.omega/LMg;
            CB.C1gag=CB.dampS*DampND;
            CB.C2gag=CB.dampNS*DampND;
            
            % ND matrices elements
            CB.M1=(CB.a-2+2/CB.a)*CB.aalpha+CB.Igag+1;
            CB.M2=CB.a*CB.aalpha+CB.Igag;
            CB.Ek=1/LMg;
            
            CB.P1 = CB.m*CB.L^2*(CB.a^2-CB.a)+CB.I;
            CB.P2 = (CB.mh+2*CB.m*(1-CB.a))*CB.L^2;
        end
        
        % %%%%%% % Get position % %%%%%% %
        function [ x, y ] = GetPos(CB, X, which)
            sS=sin(X(1)); cS=cos(X(1));
            sNS=sin(X(2)); cNS=cos(X(2));
            
            if strcmp(which,'S')
                x=CB.xS;
                y=CB.yS;
                return;
            end
            if strcmp(which,'m1')
                x=CB.xS-(1-CB.a)*CB.L*sS;
                y=CB.yS+(1-CB.a)*CB.L*cS;
                return;
            end
            if strcmp(which,'Hip')
                x=CB.xS-CB.L*sS;
                y=CB.yS+CB.L*cS;
                return;
            end
            if strcmp(which,'m2')
                x=CB.xS-CB.L*sS+CB.a*CB.L*sNS;
                y=CB.yS+CB.L*cS-CB.a*CB.L*cNS;
                return;
            end
            if strcmp(which,'NS')
                x=CB.xS-CB.L*sS+(CB.L-CB.LegShift)*sNS;
                y=CB.yS+CB.L*cS-(CB.L-CB.LegShift)*cNS;
                return;
            end
            if strcmp(which,'COM')
                % Support mass
                [ m1x, m1y ] = GetPos(CB, X, 'm1');
                % Non-support mass
                [ m2x, m2y ] = GetPos(CB, X, 'm2');
                % Hip mass
                [ mhx, mhy ] = GetPos(CB, X, 'Hip');

                x=(CB.m*m1x+CB.m*m2x+CB.mh*mhx)/(2*CB.m+CB.mh);
                y=(CB.m*m1y+CB.m*m2y+CB.mh*mhy)/(2*CB.m+CB.mh);
            end
        end
        
        % %%%%%% % Get velocity % %%%%%% %
        function [ xdot, ydot ] = GetVel(CB, X, which)
            sS=sin(X(1)); cS=cos(X(1));
            sNS=sin(X(2)); cNS=cos(X(2));
            if strcmp(which,'S')
                xdot=0;
                ydot=0;
                return;
            end
            if strcmp(which,'m1')
                xdot=-(1-CB.a)*CB.L*cS*X(3);
                ydot=-(1-CB.a)*CB.L*sS*X(3);
                return;
            end
            if strcmp(which,'Hip')
                xdot=-CB.L*cS*X(3);
                ydot=-CB.L*sS*X(3);
                return;
            end
            if strcmp(which,'m2')
                xdot=-CB.L*cS*X(3)+CB.a*CB.L*cNS*X(4);
                ydot=-CB.L*sS*X(3)+CB.a*CB.L*sNS*X(4);
                return;
            end
            if strcmp(which,'NS')
                xdot=-CB.L*cS*X(3)+(CB.L-CB.LegShift)*cNS*X(4);
                ydot=-CB.L*sS*X(3)+(CB.L-CB.LegShift)*sNS*X(4);
                return;
            end
            if strcmp(which,'COM')
                % Support mass
                [ m1x, m1y ] = GetVel(CB, X, 'm1');
                % Non-support mass
                [ m2x, m2y ] = GetVel(CB, X, 'm2');
                % Hip mass
                [ mhx, mhy ] = GetVel(CB, X, 'Hip');

                xdot=(CB.m*m1x+CB.m*m2x+CB.mh*mhx)/(2*CB.m+CB.mh);
                ydot=(CB.m*m1y+CB.m*m2y+CB.mh*mhy)/(2*CB.m+CB.mh);
            end
        end
        
        % %%%%%% % Get ground reaction forces % %%%%%% %
        function [F] = GetGRF(CB, X)
            % Calculate GRF based on the 4 DOF model (x,y,q1,q2)
            % H(q)*qdot2+h(q,qdot)*qdot=[T_stance, T_swing, Fx, Fy]T
            
            % Get qdot2
            [Xdot] = Derivative(CB, 0, X);
            q1=X(1); q2=X(2);
            q1dot=X(3); q2dot=X(4);
            q1dot2=Xdot(3); q2dot2=Xdot(4);
            
            % xdot, ydot, xdot2 and ydot2 will be set to 0 since the
            % support foot is assumed to be stationary. If the GRF
            % becomes negative then this assumption is invalid
            
            F=zeros(2,1);
            Coef1=(CB.a*CB.m-2*CB.m-CB.mh)*CB.L;
            Coef2=CB.a*CB.m*CB.L;
            F(1)=-Coef1*sin(q1)*q1dot^2 - Coef2*sin(q2)*q2dot^2 + ...
                  Coef1*cos(q1)*q1dot2 + Coef2*cos(q2)*q2dot2;
            F(2)= Coef1*cos(q1)*q1dot^2 + Coef2*cos(q2)*q2dot^2 + ...
                  Coef1*sin(q1)*q1dot2 + Coef2*sin(q2)*q2dot2 + (2*CB.m+CB.mh)*CB.g;
        end
        
        function [Length] = GetStepLength(CB, X)
            [X0,Y0]=CB.GetPos(X,'S');
            [X1,Y1]=CB.GetPos(X,'NS');
            
            Length=sqrt((X1-X0)^2+(Y1-Y0)^2);
        end
        
        function [Weight] = GetWeight(CB)
            Weight=(2*CB.m+CB.mh)*CB.g;
        end
        
        function [KE] = GetKineticEnergy(CB, X)
            [m1x, m1y]=CB.GetVel(X, 'm1');
            [mhx, mhy]=CB.GetVel(X, 'Hip');
            [m2x, m2y]=CB.GetVel(X, 'm2');
            KE = 1/2*CB.m*(m1x^2+m1y^2)+1/2*CB.mh*(mhx^2+mhy^2)+1/2*CB.m*(m2x^2+m2y^2)+1/2*CB.I*(X(3)^2+X(4)^2);
        end
        
        function [PE] = GetPotentialEnergy(CB, X)
            [m1x, m1y]=CB.GetPos(X, 'm1'); %#ok<ASGLU>
            [mhx, mhy]=CB.GetPos(X, 'Hip'); %#ok<ASGLU>
            [m2x, m2y]=CB.GetPos(X, 'm2'); %#ok<ASGLU>
            PE = CB.g*(CB.m*(m1y + m2y) + CB.mh*mhy);
        end
        
        function [Minv, N, G, Eu] = GetSystemMat(CB, X)
            sS=sin(X(1)); sNS=sin(X(2));
            sSNS=sin(X(1)-X(2)); cSNS=cos(X(1)-X(2));
            theta1t=X(3)/CB.omega;
            theta2t=X(4)/CB.omega;

            % M=[      CB.M1       -CB.aalpha*cSNS ;
            %     -CB.aalpha*cSNS        CB.M2  ];
            
            Mdet=CB.M1*CB.M2-CB.aalpha^2*cSNS^2;
            Minv=[       CB.M2       CB.aalpha*cSNS ;
                     CB.aalpha*cSNS       CB.M1  ];
            Minv=Minv/Mdet;

            N=[CB.C1gag+CB.C2gag, CB.C2gag-CB.aalpha*sSNS*theta2t;
               -CB.C2gag+CB.aalpha*sSNS*theta1t, CB.C2gag];

            G=[-(2*CB.alpha-CB.aalpha+1)*sS;
                       CB.aalpha*sNS       ];

            Eu=CB.Ek*[1 -1; 0 1]*CB.GetTorques();
        end
    
        function [U] = GetTorques(CB)
            U=[ CB.Torques(1);
                CB.Torques(2)];
        end
        
        % %%%%%% % Derivative % %%%%%% %
        function [Xdot] = Derivative(CB, t, X) %#ok<INUSL>
            % Derivative uses normal coordinates
            % Translate first to nondimensional
            Xdot=zeros(4,1);
            theta1t=X(3)/CB.omega;
            theta2t=X(4)/CB.omega;

            % Musculo-Skeletal System
            % M*quu+N*qu+G=Eu
            [Minv, N, G, Eu] = CB.GetSystemMat(X);
            quu=Minv*(Eu-N*[theta1t; theta2t]-G);
            
            % Return to normal coordinates
            qtt=CB.omega^2*quu;

            Xdot(1)=X(3); Xdot(2)=X(4);
            Xdot(3)=qtt(1); Xdot(4)=qtt(2);
        end
        
        % %%%%%% % Calculate Impact % %%%%%% %
        function [Xa] = CalcImpact(CB, Xb)
            cSNS=cos(Xb(1)-Xb(2));
            
            % Calculate Impact
            % P*q_before = Q*q_after
            P=[      CB.P1           0    ;
                CB.P2*cSNS+CB.P1   CB.P1  ];

            Q1 = -CB.a*CB.L^2*CB.m*cSNS;
            M1d = CB.M1*CB.L^2*CB.mh;
            M2d = CB.M2*CB.L^2*CB.mh;
            Q=[   Q1     M2d   ;
                Q1+M1d  Q1+M2d ];

            Xdotnew=pinv(Q)*P*[Xb(3); Xb(4)];
            % Coordinates are switched at impact,
            % angular velocities are already switched when
            % impact is calculated
            Xa=[Xb(2), Xb(1), Xdotnew(1), Xdotnew(2)];
        end
        
        % %%%%%% % Events % %%%%%% %
        function [value, isterminal, direction] = Events(CB, X, Floor)
            value=ones(CB.nEvents,1);
            isterminal=ones(CB.nEvents,1);
            direction=-ones(CB.nEvents,1);
            
            [xNS, yNS]=CB.GetPos(X,'NS');
        
            % Check for foot contact
            value(1)=yNS-Floor.Surf(xNS);
            
            % Check for robot "falling"
            [HipPosx,HipPosy]=CB.GetPos(X,'Hip'); %#ok<ASGLU>
            value(2)=HipPosy-0.7*CB.L;
        end
        
        % %%%%%% % Events % %%%%%% %
        function [CB,Xa] = HandleEvent(CB, evID, Xb, t)
            switch evID
                case 1
                    % Leg hit the ground
                    % Calculate Impact
                    Xa=CB.CalcImpact(Xb);

                    % Update support foot
                    if CB.Support==CB.Left
                        CB.Support=CB.Right;
                    else
                        if CB.Support==CB.Right
                        CB.Support=CB.Left;
                        end
                    end
                    
                    % Update support foot position
                    [xNS, yNS]=CB.GetPos(Xb,'NS');
                    CB.yS=yNS;
                    CB.xS=xNS;
        
                    % Compute velocity
                    dT = t-CB.last_t;
                    CB.curSpeed = CB.GetStepLength(Xb)/dT;
                    CB.last_t = t;
                case 2
                    % Robot fell, do nothing
                    Xa = Xb;
            end
        end
        
        function CB = SetTorques(CB,T)
            CB.Torques=T;
        end
        
        % %%%%%% % Render Compass Biped % %%%%%% %
        % In Render.m
        % Includes also DrawCircle, DrawLink and DrawVector
        
    end % end methods
end