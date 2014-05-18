classdef Simulation < handle & matlab.mixin.Copyable
    % Version 0.2 - 10/05/2014
    % This simulation integrates a system over time until
    % an event occurs, then it performs some calculations
    % and continues integrating until the next event or
    % it runs out of time
    
    properties
        Mod; % Model
        Con; % Controller
        Env; % Environment
    
        % State params
        stDim; ModCo; ConCo;
        % Event params
        nEvents; ModEv; ConEv;
    
        % Simulation parameters
        IC;
        infTime;
        tstep; tstep_normal; tstep_small = [];
        tstart; tend; tspan;
        
        % Performance tracking / Statistics
        Out; % output holder
        EndCond = 0;
        % Set EndCond to run the sim until:
        % 0 - the end of time
        % [1,numsteps] - numsteps are taken on end_slope
        % 2 - the system converges to a limit cycle
        
        CurSpeed; StepsTaken; Steps2Slope;
        MaxSlope; MinSlope;
        ICstore; nICsStored = 10;
        minDiff = 8e-7; % Min. difference for LC convergence
        stepsReq = 15; % Steps of minDiff required for convergence
        stepsSS; % Steps taken since minDiff
        
        % Poincare map calculation parameters
        IClimCyc; Period;
        PMeps = 1e-4; PMFull = 0;
        PMeigs; PMeigVs;
                
        % Rendering params
        Graphics = 1;
        Fig = 0; Once; StopSim;
        FigWidth; FigHeight; AR;
        % Environment display
        FlMin; FlMax; HeightMin; HeightMax;
        
        Follow = 1; % Follow model with camera
        % COM transformation
        tCOM; COMx0; COMy0;
        % Time display
        hTime; TimeStr = ['t = %.2f s\nOsc.=%.3f\n',...
                         'Slope = %.2f ',char(176)','\nSpeed = %s'];
        % Torques display
        nOuts; nTsteps;
        Ttime; Thold; Tbase; Tscale;
        hTorques;
        Colors = {[1 0 0],[0 0 1],[0 1 0],[0 0 0]};
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function sim = Simulation(varargin)
            switch nargin
                case 3
                    sim.Mod = varargin{1};
                    sim.Con = varargin{2};
                    sim.Env = varargin{3};
                otherwise
                    sim.Mod = CompassBiped();
                    sim.Con = Controller();
                    sim.Env = Terrain();
            end            
        end
        
        function sim = SetEndCond(sim, value)
            L = length(value);
            if L<1
                error('Invalid input for EndCond');
            end
            
            if value(1) == 1
                Error = ['When setting EndCond to 1,',...
                           'a value for num. steps is also needed',...
                           '\nPlease use sim.EndCond = [1,nsteps]'];
                if L<2
                    error(Error);
                else
                    if ~isnumeric(value(2)) || value(2)<1
                        error(Error);
                    end
                end
            end
            
            sim.EndCond = value;
        end
        
        function sim = SetTime(sim,tstart,tstep,tend)
            if nargin~=4
                error(['Set time expects 3 input arguments',...
                    ' but was provided with ',num2str(nargin)]);
            end
            sim.tstart = tstart;
            sim.tstep_normal = tstep;
            sim.tstep = tstep;
            if isnumeric(tend)
                if tend<=tstart+tstep
                    error('tend is too close to tstart');
                else
                    sim.tend = tend;
                end
                sim.infTime = 0;
            else
                if strcmp(tend,'inf')
                    % Simulation will run for indefinite time
                    sim.infTime = 1;
                    sim.tend = 10;
                end
            end
        end
        
        function [Xt] = Derivative(sim,t,X)
            sim.Mod.Torques=sim.Con.NeurOutput();

            Xt = [sim.Mod.Derivative(t,X(sim.ModCo));
                  sim.Con.Derivative(t,X(sim.ConCo))];
        end

        function [value, isterminal, direction] = Events(sim, t, X) %#ok<INUSL>
            value = zeros(sim.nEvents,1);
            isterminal = ones(sim.nEvents,1);
            direction = zeros(sim.nEvents,1);

            % Call model event function
            [value(sim.ModEv), isterminal(sim.ModEv), direction(sim.ModEv)] = ...
                sim.Mod.Events(X(sim.ModCo), sim.Env);
            % Call controller event function
            [value(sim.ConEv), isterminal(sim.ConEv), direction(sim.ConEv)] = ...
                sim.Con.Events(X(sim.ConCo));
        end

        function StopButtonCb(sim, hObject, eventdata, handles) %#ok<INUSD>
            if sim.StopSim == 0
                sim.StopSim = 1;
                sim.Out.Type = -1;
                sim.Out.Text = 'Simulation stopped by user';
                set(hObject,'String','Close Window');
            else
                close(sim.Fig)
            end
        end  % StopButtonCallback
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ControlParams=[1.266646286756660, 0.597267606409947,...
%                -7.384173104761088, 0.126800300703305, 0.072267297914989,...
%                 5.191278913625541, 0.166498744328709, 0.053696788954836];
%                     
% SimType=1;
% 
% if SimType==1
%     % Terrain parameters
%     FloorType=0;
%     InitSlope=0;
%     EndSlope=0;
%     Adaptive=0;
% else
%     % Terrain parameters
%     FloorType=2;
%     InitSlope=0;
%     Adaptive=1;
%     
%     % Start the slope after giving the robot
%     % 10 seconds to reach steady state
%     Floor.start_x=AvgVel*10;
%     
%     MinSlope=0;
%     MaxSlope=0;
% end
% 
% 
%     function StopButtonCallback(t,X) %#ok<INUSD>
%         StopSim=1;
%         close all
%         dispLength='%5.7f';
%         disp(['InitCond=[',num2str(IC_Store(1,1),dispLength),', ',num2str(IC_Store(2,1),dispLength),', ',...
%                         num2str(IC_Store(3,1),dispLength),', ',num2str(IC_Store(4,1),dispLength),', ',...
%                         num2str(IC_Store(5,1),dispLength),'];',10]);
%     end  % StopButtonCallback

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

