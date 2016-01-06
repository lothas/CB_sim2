function MatsuokaTest()
%MATSUOKATEST Summary of this function goes here
%   Detailed explanation goes here

MO = Matsuoka();
MO.tau = 0.3;
MO.tav = 1;     % Affects oscillations period and shape
MO.beta = 3;   % Affects oscillations period (must be close to wfe?)
MO.u0 = 3;        % Affects final amplitude
MO.wfe = -3;

uSS = MO.u0/(1+MO.beta);

% Initial conditions
IC0 = [uSS; uSS; -uSS; 0];

% tspan
t_start = 0;
t_end = 10;
t_step = 0.05;
t_span = t_start:t_step:t_end;

% Run simulation
options = odeset('MaxStep',t_step/10,'RelTol',.5e-7,'AbsTol',.5e-8);
[TTemp,XTemp] = ode45(@MO.Derivative,t_span,IC0,options);
% options = odeset('MaxStep',t_step/10,'RelTol',.5e-7,'AbsTol',.5e-8,...
%             'Events', @MO.Events);
% [TTemp,XTemp,TE,YE,IE] = ...
%     ode45(@MO.Derivative,t_span,IC0,options); %#ok<ASGLU>

y = max(XTemp(:,1),0) - 0.1*max(XTemp(:,3),0);
subplot(2,1,1)
plot(TTemp, XTemp(:,[1,3]));
subplot(2,1,2)
plot(TTemp, y);
grid on

end

