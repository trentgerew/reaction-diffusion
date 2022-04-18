clear
clc
close all

%% Parameters
load('params')

%% Simulation
% Form matrix
M = [-k12 , k21 ; k12 , -k21];

% Solve
tspan = [0 tstop]; % time domain
[t,p] = ode45(@(t,p) odefun(t,p,M), tspan, init); % solve linear system
p = p'; % tranpose solution so rows are n and cols are t

% Plot
figure;
plot(t,p(1,:),'r-',t,p(2,:),'b-')
xlabel('t')
legend('\rho_1(t)','\rho_2(t)')

% Show end probability
formatSpec = 'rho_%d (t_f) = %1.8f\n';
fprintf(formatSpec,1,p(1,end))
fprintf(formatSpec,2,p(2,end))

%% Functions
function dpdt = odefun(~,p,M)
    dpdt = M * p;
end