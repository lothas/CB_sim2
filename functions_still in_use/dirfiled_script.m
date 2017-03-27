% dirfiled.m script
close all; clear all; clc;

%% % define the ODE
% Parameters ranges:
t_min = 0;
dt = 0.5;
t_max = 10;
x_min = 0;
dx = 0.5;
x_max = 10;

syms w x t
f00(x,t) = (tanh(w*x)-x);
w0 = 0.1:0.1:1.5;


for i = 1:length(w0)
    f0=subs(f00,w,w0(1,i));
    f=matlabFunction(f0);
    dirfield(f,...
        t_min:dt:t_max,...
        x_min:dx:x_max);
    pause(0.1);
    cla
end

%%
close all; clear all; clc;
I = 0.1;
u0 = 1;
m = 1;
g = 9.81;
l = 2;
b = 0.01;
syms Y(t)
[V] = odeToVectorField(I*diff(Y, 2) == u0 - m*g*l*sin(Y) - b*diff(Y));
M = matlabFunction(V,'vars', {'t','Y'});
sol = ode45(M,[0 20],[2 0]);

time = sol.x;
theta = sol.y(1,:);
theta_dot = sol.y(2,:);
figure;
plot(theta,theta_dot);