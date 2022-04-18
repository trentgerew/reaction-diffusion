clear
clc
close all

%% Simulation Parameters
tstop = .001; % stopping time

%% Problem Parameters
% Rate constants
k12 = 2000;
k21 = 1000;

% Initial conditions
p0 = [1 0];

%% Simulation
% Form matrix
M = [-k12 , k21 ; k12 , -k21];

% Solve
tspan = [0 tstop]; % time domain
[t,p] = ode45(@(t,p) odefun(t,p,M), tspan, p0); % solve linear system
p = p'; % tranpose solution so rows are n and cols are t

% Plot
figure;
plot(t,p(1,:),'r-',t,p(2,:),'b-')
xlabel('t')
legend('\rho_1(t)','\rho_2(t)')

%% Functions
function dpdt = odefun(~,p,M)
    dpdt = M * p;
end