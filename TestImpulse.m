function TestImpulse()
Graphics = 0;

Mod = CompassBiped();
Mod = Mod.Set('damp',0,'I',0); % no damping, point-mass

% Gait parameters
alpha = 0.1126761;
theta_dot = [-0.5209123, -0.5055816];
delta_theta_dot = [-0.02,  -0.05];

% Different pulse durations to try
dts = [0.001, 0.003, 0.01, 0.03, 0.1];
    
filename = 'Con_Analysis.mat';
if exist(filename,'file') ~= 2
    % Run simulation with impulsive controller
    Sim1 = ImpulseSim(alpha, theta_dot, delta_theta_dot);

    res(length(dts)+1) = AnalyzeContr(Sim1);

    % Do some plots
    disp(res(end).eig_val);

    frmt_length = 7;
    disp(['alpha = ',num2str(res(end).Sim.ICstore(1,1), frmt_length),';',10,...
        'theta_dot = [',num2str(res(end).Sim.ICstore(3,1), frmt_length),', '...
                        num2str(res(end).Sim.ICstore(4,1), frmt_length),'];']) %;,10,...
    %     'delta = [',num2str(Sim1.ICstore(3,1)),', '...
    %                     num2str(Sim1.ICstore(4,1)),'];']);

    % Model data
    m = Mod.m; mh = Mod.mh; L = Mod.L; IL = Mod.I; a = Mod.a;
    M = zeros(2,2);
    M(1,1) = IL + L^2*(m*a^2 - 2*m*a + mh + 2*m);
    M(1,2) = -L^2*a*m*cos(2*alpha);
    M(2,1) = M(1,2);
    M(2,2) = m*L^2*a^2 + IL;
    E = [1 -1; 0 1];

    % Calculate required impulse to obtain the difference in ang. vel.
    delta = (E^-1*M*delta_theta_dot')';

    for i = 1:length(dts)
        % Run simulations with quasi-impulsive controller with
        % different dt
        Sim = ImpulseSim(alpha, theta_dot, delta, dts(i));
        try
            res(i) = AnalyzeContr(Sim);
        catch
            disp(['Failed to analyze dt = ',num2str(dts(i))]);
        end
    end

    save(filename,'delta_theta_dot','res');
else
    % Plot results
    data = load(filename);
    legends = cell(length(data.res),1);
    h = zeros(length(data.res),1);
    close all
    
    figure
    hold on
    for i = 1:length(data.res)-1
        if ~isempty(data.res(i).eig_val)
            h(i) = PlotRes(data.res(i));
        end
        legends{i} = num2str(dts(i));
    end
    h(end) = PlotRes(data.res(end),'-',2);
    legends{end} = 'impulsive';
    
    % Clear handle+legend of failed simulations
    out = find(h==0);
    h(out) = [];
    legends(out) = [];
    
    % Add legend
    legend(h, legends);
end

    function res = AnalyzeContr(Sim, do_plot)
        if nargin<3
            do_plot = false;
        end
        
        Sim.Graphics = Graphics;
        aSim = deepcopy(Sim); % save clean copy

        % Simulate
        Sim = Sim.Run();
        
        if Sim.Out.Type == 5 % Sim converged, analyze
            % Calculate eigenvalues
            [res.eig_val, ~] = Sim.Poincare();
        
            % Get limit cycle and events
            aSim.IC = Sim.IClimCyc;
            aSim.SetTime(0,0.005,2); % smaller time step
            aSim.EndCond = [1,1]; % stop after 1 step
            aSim.Graphics = Graphics;
            aSim = aSim.Run();

            % Get state right before impact, right after impact and
            % after adding the impulse
            Xminus = aSim.Out.X(end,:);
            Xplus_imp = aSim.Out.X(1,:);
            Xplus = [aSim.Mod.CalcImpact(Xminus), Xminus(aSim.ConCo)];
            
            % Prepare output structure
            res.Sim = aSim;
            res.delta = Xplus_imp(3:4) - Xplus(3:4);
            res.alpha = (aSim.Out.X(1,1) - aSim.Out.X(1,2))/2;
            res.period = aSim.Out.T(end);
            res.X = aSim.Out.X;
            res.T = aSim.Out.T;
            res.discnt = [Xminus; Xplus; Xplus_imp];

            if do_plot
                figure
                hold on
                PlotRes(res);
            end
        end
    end

    function h = PlotRes(res, style, lw)
        if nargin == 1
            style = '--';
            lw = 1;
        end
        h = plot([res.X(:,1);res.X(:,2)],[res.X(:,3);res.X(:,4)], ...
            style,'LineWidth',lw);
        scatter(res.discnt(1, [1,2]), res.discnt(1, [3,4]), 'o')
        scatter(res.discnt(2, [1,2]), res.discnt(2, [3,4]), '+')
        scatter(res.discnt(3, [1,2]), res.discnt(3, [3,4]), '^')
        title(['Limit cycle for:  dth_1_t = ', ...
                num2str(res.delta(1), 4), ...
                '  |  dth_2_t = ', ...
                num2str(res.delta(2), 4)]);
        xlabel(['alpha = ', num2str(res.alpha, 4), ...
                '  |  T = ', num2str(res.period, 4)]);
    end
end