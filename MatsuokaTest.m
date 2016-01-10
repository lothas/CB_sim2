function MatsuokaTest()
%MATSUOKATEST Summary of this function goes here
%   Detailed explanation goes here

MO = Matsuoka();
MO.nPulses = 2;
MO.tau = 0.5;
MO.tav = 0.25;     % Affects oscillations period and shape
MO.beta = 25;   % Affects oscillations period (must be close to wfe?)
MO.u0 = 1;        % Affects final amplitude
wfe = -5;
D = mat2cell(repmat([0, wfe; wfe, 0],1,MO.nPulses), ...
             2,2*ones(1,MO.nPulses)); 
MO.win = blkdiag(D{:});
% MO.win = blkdiag([0, -5; -25, 0],[0, -10; -10, 0]);
% MO.wex = zeros(2*MO.nPulses);
MO.wex = [ 0,  0,  0, -15;
           0,  0,  0,  0;
           0,  0,  0,  0;
           -15,  0,  0,  0];
MO.OutM = [ 1, -1,  0,  0;
            0,  0,  1, -1];

% Initial conditions
u_ss = MO.u0/MO.beta*abs(wfe);
IC0 = [u_ss, 0, u_ss, 0, 0, 0, 0, 0];

% tspan
t_start = 0;
t_end = 10;
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
figure
N = MO.nPulses;
for i = 1:4
    subplot(4,1,i)
    plot(TTemp, XTemp(:,N*(i-1)+1:N*i));
end
figure
plot(TTemp, MO.Output(0, XTemp', 0));
grid on

end

