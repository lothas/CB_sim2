function TestImpulse()
Graphics = 0;

% Gait parameters
% alpha = 0.1126761;
% theta_dot = [-0.5209123, -0.5055816];
% delta_theta_dot = [-0.02,  -0.05];

alpha = 0.05828347;
theta_dot = [-0.2718532, -0.1238651];
delta_theta_dot = 0.1*[0.0871, 1.5031];

% Different pulse durations to try
dts = [0.001, 0.003, 0.01, 0.03, 0.1];
    
filename = 'Con_Analysis.mat';
if exist(filename,'file') ~= 2
    % Run simulation with impulsive controller
    Sim1 = ImpulseSim(alpha, theta_dot, delta_theta_dot);

    res(length(dts)+1) = AnalyzeContr(Sim1);
    
    if ~isfield(res(end), 'eig_val')
        % This will fail if the simulation didn't converge

        frmt_length = 7;
        disp(['alpha = ',num2str(res(end).Sim.ICstore(1,1), frmt_length),';',10,...
            'theta_dot = [',num2str(res(end).Sim.ICstore(3,1), frmt_length),', '...
                            num2str(res(end).Sim.ICstore(4,1), frmt_length),'];']) %;,10,...
        %     'delta = [',num2str(Sim1.ICstore(3,1)),', '...
        %                     num2str(Sim1.ICstore(4,1)),'];']);
        
        return
    end
    
    % Do some plots
    disp(res(end).eig_val);

    % Model data
    m = Sim1.Mod.m; mh = Sim1.Mod.mh; L = Sim1.Mod.L;
    IL = Sim1.Mod.I; a = Sim1.Mod.a;
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
    set(gca,'FontSize',14);
    
    % Plot eigenvalues
    eig_val = [data.res.eig_val];
    eig_val = [eig_val(:,end), eig_val(:,1:end-1)];
    dts = zeros(1, size(eig_val,2));
    for i = 2:length(dts)
        dts(i) = mean(data.res(i-1).Sim.Con.Duration);
    end
    
    figure
    h1 = semilogx(dts, abs(eig_val)','*-.','LineWidth',2);
    hold on
    for i = 1:size(eig_val,1)
        h2 = plot(xlim,repmat(abs(eig_val(i,1)),1,2),'k--');
    end
%     plot(dts, abs(eig_val)','*-.','LineWidth',2);
    ylabel('|\lambda_i|','FontSize',16);
    xlabel('\Deltat','FontSize',16);
    set(gca,'FontSize',14);
    
    legend([h1(1),h2],{'Quasi-impulsive','Impulsive'})
    
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
            Xplus_imp = GetStateAfterImp(aSim);
            Xplus = [aSim.Mod.CalcImpact(Xminus), Xminus(aSim.ConCo)];
            
            % Prepare output structure
            res.Sim = aSim;
            res.delta = Xplus_imp(3:4) - Xplus(3:4);
            res.alpha = (aSim.Out.X(1,1) - aSim.Out.X(1,2))/2;
            res.period = aSim.Out.T(end);
            if aSim.Con.FBImpulse>0
                % Add the impact state to the LC
                res.X = [aSim.Out.X; Xplus([2 1 4 3 5])];
                res.T = [aSim.Out.T; aSim.Out.T(end)];
            else
                res.X = aSim.Out.X;
                res.T = aSim.Out.T;
            end
            res.discnt = [Xminus; Xplus; Xplus_imp];

            if do_plot
                figure
                hold on
                PlotRes(res);
            end
        else
            res.Sim = Sim;
        end
    end

    function x_imp = GetStateAfterImp(sim)
        if sim.Con.FBImpulse>0
            % Impulsive case
            x_imp = sim.Out.X(1,:);
        else
            if ~isempty(sim.Con.ExtPulses)
                % Quasi-impulsive case
                T_end = ind2sub(size(sim.Out.Torques), ...
                    find(sim.Out.Torques == 0, 1, 'first'));
                x_imp = sim.Out.X(T_end,:);
            else
                x_imp = sim.Out.X(1,:);
            end
        end
    end

    function h = PlotRes(res, style, lw)
        MarkerSize = 50;
        LineWidth = 1.5;
        if nargin == 1
            style = '--';
            lw = 1;
        end
        h = plot([res.X(:,1);res.X(:,2);res.X(1,1)], ...
                 [res.X(:,3);res.X(:,4);res.X(1,3)], ...
            style,'LineWidth',LineWidth+lw);
        scatter(res.discnt(1, [1,2]), res.discnt(1, [3,4]), MarkerSize, ...
            'o', 'LineWidth', LineWidth);
        scatter(res.discnt(2, [1,2]), res.discnt(2, [3,4]), MarkerSize, ...
            '+', 'LineWidth', LineWidth);
        scatter(res.discnt(3, [1,2]), res.discnt(3, [3,4]), MarkerSize, ...
            '^', 'LineWidth', LineWidth);
        title(['Limit cycle for:  dth_1_t = ', ...
                num2str(res.delta(1), 4), ...
                '  |  dth_2_t = ', ...
                num2str(res.delta(2), 4)],...
                'FontSize',16);
        xlabel(['alpha = ', num2str(res.alpha, 4), ...
                '  |  T = ', num2str(res.period, 4)],...
                'FontSize',16);
    end
end