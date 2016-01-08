function MatsuokaTest()
%MATSUOKATEST Summary of this function goes here
%   Detailed explanation goes here

MO = Matsuoka();
MO.tau = 0.25;
MO.tav = 0.5;     % Affects oscillations period and shape
MO.beta = 20;   % Affects oscillations period (must be close to wfe?)
MO.u0 = 60;        % Affects final amplitude
MO.wfe = -5;

uSS = MO.u0/(1-MO.wfe);

% Initial conditions
IC0 = [-15.7618    3.3300    2.7334    6.8556]/2;

% tspan
t_start = 0;
t_end = 0.754*2;
t_step = 0.03;
t_span = t_start:t_step:t_end;

% Run simulation
options = odeset('MaxStep',t_step/10,'RelTol',.5e-7,'AbsTol',.5e-8);
[TTemp,XTemp] = ode45(@MO.Derivative,t_span,IC0,options);
% options = odeset('MaxStep',t_step/10,'RelTol',.5e-7,'AbsTol',.5e-8,...
%             'Events', @MO.Events);
% [TTemp,XTemp,TE,YE,IE] = ...
%     ode45(@MO.Derivative,t_span,IC0,options); %#ok<ASGLU>

y = max(XTemp(:,1:2:end),0);
out = y(:,1) - 0.1*y(:,2);
figure
for i = 1:4
    subplot(4,1,i)
    plot(TTemp, XTemp(:,i));
end
figure
plot(TTemp, out);
grid on

end

