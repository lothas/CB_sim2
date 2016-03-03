run_imp = 1;

Mod = CompassBiped();
Mod = Mod.Set('damp',0,'I',0); % no damping, point-mass

% Gait parameters
alpha = 0.1126761;
theta_dot = [-0.5209123, -0.5055816];
delta_theta_dot = [-0.02,  -0.05];
    
if run_imp == 1
    Sim1 = ImpulseSim(alpha, theta_dot, delta_theta_dot);
    aSim1 = deepcopy(Sim1);
    
    % Simulate
    Sim1 = Sim1.Run();

    % Calculate eigenvalues
    [EigVal,EigVec] = Sim1.Poincare();

    % Do some plots
    disp(EigVal);
    
    frmt_length = 7;
    disp(['alpha = ',num2str(Sim1.ICstore(1,1), frmt_length),';',10,...
        'theta_dot = [',num2str(Sim1.ICstore(3,1), frmt_length),', '...
                        num2str(Sim1.ICstore(4,1), frmt_length),'];']) %;,10,...
    %     'delta = [',num2str(Sim1.ICstore(3,1)),', '...
    %                     num2str(Sim1.ICstore(4,1)),'];']);
    
    if Sim1.Out.Type == 5
        % Sim converged, analyze
        aSim1.IC = Sim1.IClimCyc;
        aSim1.SetTime(0,0.001,2);
        aSim1.EndCond = [1,1];
        aSim1 = aSim1.Run();
        
        Xminus = aSim1.Out.X(end,:);
        Xplus_imp = aSim1.Out.X(1,:);
        Xplus = [aSim1.Mod.CalcImpact(Xminus), Xminus(aSim.ConCo)];
        
        alpha = (aSim1.Out.X(1,1) - aSim1.Out.X(1,2))/2;
        T = aSim1.Out.T(end);
        
        figure()
        plot(aSim1.Out.X(:,1),aSim1.Out.X(:,3));
        hold on
        plot(aSim1.Out.X(:,2),aSim1.Out.X(:,4));
        scatter(Xminus([1,2]), Xminus([3,4]), 'o')
        scatter(Xplus_imp([1,2]), Xplus_imp([3,4]), '^')
        scatter(Xplus([1,2]), Xplus([3,4]), '+')
        title(['Limit cycle for:  dth_1_t = ', ...
                num2str(delta_theta_dot(1), 4), '  |  dth_2_t = ', ...
                num2str(delta_theta_dot(2), 4)]);
        xlabel(['alpha = ', num2str(alpha, 4), ...
                '  |  T = ', num2str(T, 4)]);
        
    end
end

% Model data
m = Mod.m; mh = Mod.mh; L = Mod.L; IL = Mod.I; a = Mod.a;
M = zeros(2,2);
M(1,1) = IL + L^2*(m*a^2 - 2*m*a + mh + 2*m);
M(1,2) = -L^2*a*m*cos(2*alpha);
M(2,1) = M(1,2);
M(2,2) = m*L^2*a^2 + IL;
E = [1 -1; 0 1];

delta = (E^-1*M*delta_theta_dot')';
dt = 0.01;

Sim2 = ImpulseSim(alpha, theta_dot, delta, dt);